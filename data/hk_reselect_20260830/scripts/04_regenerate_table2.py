#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
重算 Table 2 的 GTEx ACAT-O 列（与新看家基因同一管道），eQTLGen 列保留已发表值。
"""
import io, sys, os
import pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

ROOT = r'E:\workbuddy\hk_reselect_20260830'
PROC = r'E:\workbuddy\eqtl-source-discordance-audit\data\processed'

# 1) 已发表 Table 2 (eQTLGen 列保留)
pub = pd.read_csv(os.path.join(PROC, 'enrichment_comparison.csv'))
# 2) 同一管道的 GTEx ACAT-O 重算（含 HK_v2）
arms = pd.read_csv(os.path.join(ROOT, 'arms_all_groups.csv'))
# arms.grp 用的是 '30 HOTAIR Candidate'，pub.Group_eqtl 是 '30_HOTAIR_Candidate'，做映射
gmap = {'30 HOTAIR Candidate': '30_HOTAIR_Candidate',
        '44 Non-Candidate': '44_NonCandidate_HOTAIR',
        '30 T2DM Control': '30_T2DM_Control',
        'HK_v2': 'Housekeeping v2'}

import numpy as np, math
def bh_q(p):
    p = np.asarray(p, float); m = len(p)
    if m == 0: return p
    order = np.argsort(p)
    qs = p[order] * m / np.arange(1, m + 1)
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    q = np.empty(m); q[order] = np.minimum(qs, 1.0); return q

new = []
for (g, t), sub in arms.groupby(['grp', 'trait']):
    s = sub.dropna(subset=['p_ACAT_O'])
    n = len(s)
    q = bh_q(s['p_ACAT_O'].values)
    n_fdr = int((q < 0.05).sum())
    new.append({'Group_internal': gmap.get(g, g), 'Trait': t,
                'N_Tested_GTEx_new': n, 'N_FDR_GTEx_new': n_fdr,
                'GTEx_FDR_new': round(100 * n_fdr / n, 1)})
new = pd.DataFrame(new)

# 合并
mrg = pub.merge(new, left_on=['Group_eqtl', 'Trait'],
                right_on=['Group_internal', 'Trait'], how='left')

# 拼装新的 Table 2（GTEx ACAT-O 列替换，eQTLGen 列保留）
out = mrg[['Trait', 'Group_gtex']].copy()
out = out.rename(columns={'Group_gtex': 'Gene Group'})
out['N Tested (GTEx)'] = mrg['N_Tested_GTEx_new'].astype('Int64')
out['N FDR (GTEx)'] = mrg['N_FDR_GTEx_new'].astype('Int64')
out['GTEx ACAT-O FDR (%)'] = mrg['GTEx_FDR_new']
out['N Tested (eQTLGen)'] = mrg['N_Tested_eqtl'].astype('Int64')
out['N FDR (eQTLGen)'] = mrg['N_FDR005_eqtl'].astype('Int64')
out['eQTLGen BH FDR (%)'] = mrg['eQTLGen_FDR%']

# 排序：与已发表一致
order = [('30_HOTAIR_Candidate', 'DR'), ('30_HOTAIR_Candidate', 'DN'), ('30_HOTAIR_Candidate', 'DPN'),
         ('44_NonCandidate_HOTAIR', 'DR'), ('44_NonCandidate_HOTAIR', 'DN'), ('44_NonCandidate_HOTAIR', 'DPN'),
         ('30_T2DM_Control', 'DR'), ('30_T2DM_Control', 'DN'), ('30_T2DM_Control', 'DPN')]
out['__k'] = out.apply(lambda r: order.index((r['Gene Group'], r['Trait'])), axis=1)
out = out.sort_values('__k').drop(columns='__k').reset_index(drop=True)

# 添加 Housekeeping v2
hk = new[new.Group_internal == 'Housekeeping v2'].copy()
hk = hk.assign(**{'Gene Group': 'Housekeeping v2',
                  'N Tested (GTEx)': hk['N_Tested_GTEx_new'].astype('Int64'),
                  'N FDR (GTEx)': hk['N_FDR_GTEx_new'].astype('Int64'),
                  'GTEx ACAT-O FDR (%)': hk['GTEx_FDR_new'],
                  'N Tested (eQTLGen)': pd.NA,
                  'N FDR (eQTLGen)': pd.NA,
                  'eQTLGen BH FDR (%)': pd.NA})[list(out.columns)]
out = pd.concat([out, hk], ignore_index=True)

print('新 Table 2（GTEx ACAT-O 列由本次同一管道重算，eQTLGen 列保留已发表）:')
print(out.to_string(index=False))
out.to_csv(os.path.join(ROOT, 'Table2_v2.csv'), index=False)

# 旧/新对照
cmp = mrg[['Group_gtex', 'Trait',
           'GTEx_FDR%', 'GTEx_FDR_new']].copy()
cmp = cmp.rename(columns={'Group_gtex': 'Group', 'GTEx_FDR%': 'Published %', 'GTEx_FDR_new': 'Re-run %'})
cmp['delta_pp'] = (cmp['Re-run %'] - cmp['Published %']).round(1)
print('\n旧/新 GTEx ACAT-O 对照:')
print(cmp.to_string(index=False))
cmp.to_csv(os.path.join(ROOT, 'Table2_published_vs_rerun.csv'), index=False)
