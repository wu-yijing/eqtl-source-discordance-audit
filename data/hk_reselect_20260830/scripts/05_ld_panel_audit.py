#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LD 面板审计：定量评估 g1000_eur（22,665,064 SNPs）覆盖不足导致的基因损失，
并核对已发表 per-tissue CSV 中的不可复现行。
"""
import io, os, sqlite3, sys, csv, math
import pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

PROC = r'E:\workbuddy\eqtl-source-discordance-audit\data\processed'
MODEL_DIR = r'E:\workbuddy\BMC Genomics投稿资料\DR，DN，DM芬兰原始数据\mashr_eqtl\eqtl\mashr'
BED_PREFIX = r'E:\workbuddy\2026-06-24-05-57-20\tools\ldref\g1000_eur'
TISSUES = ['Nerve_Tibial', 'Whole_Blood']
PHENOS = ['DR', 'DN', 'DPN']

cov = pd.read_csv(os.path.join(PROC, 'covariate_matrix.csv'))
cov['G'] = cov.Gene.str.upper()
grp = dict(zip(cov['G'], cov['Group']))
panel = sorted(grp)

# 1) 收集所有模型 SNP 与"面板缺失"情况
model_snps = {}          # (tissue, gene) -> [(rsid, weight)]
for t in TISSUES:
    conn = sqlite3.connect(os.path.join(MODEL_DIR, f'mashr_{t}.db'))
    meta = {}
    for g, gn, n in conn.execute('SELECT gene, genename, "n.snps.in.model" FROM extra'):
        if gn:
            meta[gn.upper()] = g
    for gname, gid in meta.items():
        if gname not in grp:
            continue
        w = conn.execute('SELECT rsid, weight FROM weights WHERE gene=?', (gid,)).fetchall()
        if w:
            model_snps[(t, gname)] = w
    conn.close()

wanted = {rs for w in model_snps.values() for rs, _ in w}
bed_ids = set()
with open(BED_PREFIX + '.bim') as f:
    for line in f:
        bed_ids.add(line.split('\t', 2)[1])
missing = wanted - bed_ids
print(f'104 panel 模型 SNP 总数（去重）: {len(wanted)}')
print(f'  LD 面板中缺失: {len(missing)}  ->  {sorted(missing)}')

# 2) 逐基因判断：是否因缺失 SNP 而无法计算
lost = []
for (t, g), w in model_snps.items():
    bad = [rs for rs, _ in w if rs not in bed_ids]
    if bad:
        lost.append({'tissue': t, 'gene': g, 'group': grp[g],
                     'n_model': len(w), 'n_missing': len(bad), 'missing': ', '.join(bad)})
lost = pd.DataFrame(lost).sort_values(['group', 'gene', 'tissue'])
print(f'\n因 LD 面板缺 SNP 而无法计算的 gene×tissue 组合: {len(lost)}')
print(lost.to_string(index=False))

by_grp = lost.groupby(['group', 'tissue']).size().unstack(fill_value=0)
print('\n按分组统计可计算性损失（gene×tissue）:')
print(by_grp.to_string())

# 3) 与已发表 per-tissue CSV 核对
print('\n' + '=' * 78)
print('已发表 per-tissue CSV 的可复现性核对')
print('=' * 78)
mine = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'validate_panel104.csv'))
for t in TISSUES:
    for ph in PHENOS:
        ref = pd.read_csv(os.path.join(PROC, f'gtex_{t}_{ph}.csv'))
        m = mine[(mine.tissue == t) & (mine.trait == ph)]
        a, b = set(ref.gene.str.upper()), set(m.gene.str.upper())
        extra = sorted(a - b)
        extra_panel = [g for g in extra if g in grp]
        extra_other = [g for g in extra if g not in grp]
        print(f'{t:<13} {ph:<4} 已发表 {len(a):>3}  可复现 {len(b):>3}  '
              f'不可复现 {len(extra):>2} (panel内 {len(extra_panel)}, panel外 {len(extra_other)})')
        if ph == 'DR':
            print(f'   panel内不可复现: {extra_panel}')
            print(f'   panel外不可复现: {extra_other}')

# 4) 各分组可测基因数：已发表 ACAT-O 文件 vs 本次重跑
print('\n' + '=' * 78)
print('各分组 GTEx ACAT-O 可测基因数三方对照')
print('=' * 78)
aco = pd.read_csv(os.path.join(PROC, 'gtex_acat_o_results.csv'))
aco['G'] = aco.Gene.str.upper()
aco['grp'] = aco.G.map(grp)
arms = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'arms_all_groups.csv'))
rows = []
for g in ['30 HOTAIR Candidate', '44 Non-Candidate', '30 T2DM Control']:
    n_aco = int((aco[aco.grp == g].groupby('Trait').size()).mean())
    n_new = int((arms[arms.grp == g].dropna(subset=['p_ACAT_O']).groupby('trait').size()).mean())
    n_total = sum(1 for x in grp.values() if x == g)
    rows.append({'Group': g, '组内基因数': n_total,
                 '已发表ACAT-O文件': n_aco, '本次重跑': n_new})
print(pd.DataFrame(rows).to_string(index=False))
print('\n注：已发表 Table 2 的 "N Tested (GTEx)" 列为 27 / 33 / 19，')
print('    既不等于归档 ACAT-O 文件（27 / 27 / 19），也不等于本次重跑（26 / 26 / 19）。')
