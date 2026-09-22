#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HK v2 全量落盘校验（verify_all.py）

四重检查：
  A. 落盘清单   —— 所有产出是否存在、大小、修改时间
  B. 文件完整性 —— docx 是否为有效 OOXML（zip 无损、可打开、段落/表格/图片数）
  C. 数字一致性 —— 关键统计量在 hk_v2_summary.json / 英文版 / 中文版 / 两份 MD 之间逐一对齐
  D. 无残留     —— v1 污染披露与"前稿"表述在两份手稿中零命中

只读脚本，不修改任何文件。
"""
import io, os, sys, json, re, zipfile, datetime

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

BASE = r'E:/workbuddy/BMC Genomics投稿资料'
ROOT = r'E:/workbuddy/hk_reselect_20260830'

EN      = os.path.join(BASE, '精修手稿19_HKv2.docx')
ZH      = os.path.join(BASE, '中文精修手稿19_HKv2.docx')
MD_SUM  = os.path.join(ROOT, 'HKv2_交付汇总.md')
MD_DIFF = os.path.join(ROOT, '结论差异分析_HKv1vsV2.md')
GENES   = os.path.join(ROOT, 'hk_genes_v2.txt')
SUMMARY = os.path.join(ROOT, 'hk_v2_summary.json')
SCRIPT  = os.path.join(ROOT, 'zh_translate', 'strip_v1_disclosure.py')

RESULTS = []


def rec(section, item, ok, note=''):
    RESULTS.append((section, item, ok, note))
    print(f"  [{'PASS' if ok else 'FAIL'}] {item}" + (f"  — {note}" if note else ''))


def ts(p):
    return datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime('%Y-%m-%d %H:%M:%S')


# ============================================================ A. 落盘清单
print('\n' + '=' * 78)
print('A. 落盘清单（存在性 / 大小 / 修改时间）')
print('=' * 78)
FILES = [
    ('英文版手稿', EN, True),
    ('中文版手稿', ZH, True),
    ('交付汇总 MD', MD_SUM, True),
    ('结论差异分析 MD', MD_DIFF, True),
    ('看家基因权威清单', GENES, True),
    ('统计汇总 JSON', SUMMARY, True),
    ('去披露脚本', SCRIPT, True),
    ('重跑: arms_all_groups.csv', os.path.join(ROOT, 'arms_all_groups.csv'), False),
    ('重跑: TableS6_hk_control_v2.csv', os.path.join(ROOT, 'TableS6_hk_control_v2.csv'), False),
    ('重跑: Figure9_hk_v2.png', os.path.join(ROOT, 'Figure9_hk_v2.png'), False),
    ('重跑: Figure9_hk_v2.pdf', os.path.join(ROOT, 'Figure9_hk_v2.pdf'), False),
]
for label, path, required in FILES:
    if os.path.exists(path):
        size = os.path.getsize(path)
        print(f"  [PASS] {label:<26} {size:>9,} B   {ts(path)}")
        rec('A', label, True, f'{size:,} B')
    else:
        print(f"  [FAIL] {label:<26} 缺失: {path}")
        rec('A', label, not required, '缺失')

# 备份文件
print('\n  备份文件：')
baks = [f for f in os.listdir(BASE) if f.startswith('_backup_') and 'StripV1' in f]
for b in sorted(baks):
    print(f"    · {b}  ({os.path.getsize(os.path.join(BASE,b)):,} B)")
rec('A', '修改前备份（2 份）', len(baks) >= 2, f'找到 {len(baks)} 份')


# ============================================================ B. 文件完整性
print('\n' + '=' * 78)
print('B. 文件完整性（OOXML 有效性 / 结构计数）')
print('=' * 78)
from docx import Document

def docx_integrity(path, label):
    # zip 层
    try:
        z = zipfile.ZipFile(path)
        bad = z.testzip()
        z.close()
        if bad is None:
            rec('B', f'{label} zip 完整性', True, '无损坏条目')
        else:
            rec('B', f'{label} zip 完整性', False, f'损坏: {bad}')
            return None
    except Exception as e:
        rec('B', f'{label} zip 完整性', False, str(e))
        return None
    # python-docx 层
    try:
        d = Document(path)
        n_par = len(d.paragraphs)
        n_tab = len(d.tables)
        n_img = sum(1 for r in d.part.rels.values() if 'image' in r.reltype)
        expect = (n_par == 254 and n_tab == 3 and n_img == 10)
        rec('B', f'{label} 结构 254段/3表/10图', expect, f'实得 {n_par}段/{n_tab}表/{n_img}图')
        return d
    except Exception as e:
        rec('B', f'{label} 可被 python-docx 打开', False, str(e))
        return None

doc_en = docx_integrity(EN, '英文版')
doc_zh = docx_integrity(ZH, '中文版')


# ============================================================ C. 数字一致性
print('\n' + '=' * 78)
print('C. 数字一致性（JSON / 英文版 / 中文版 / 两份 MD 五方对齐）')
print('=' * 78)

with open(SUMMARY, encoding='utf-8') as f:
    summ = json.load(f)
with open(GENES, encoding='utf-8') as f:
    truth_genes = [l.strip() for l in f if l.strip() and not l.startswith('#')]

en_text = '\n'.join(p.text for p in doc_en.paragraphs) if doc_en else ''
zh_text = '\n'.join(p.text for p in doc_zh.paragraphs) if doc_zh else ''
md_sum  = open(MD_SUM,  encoding='utf-8').read() if os.path.exists(MD_SUM)  else ''
md_diff = open(MD_DIFF, encoding='utf-8').read() if os.path.exists(MD_DIFF) else ''

# C1 基因集
print('\n  C1. 看家基因集（30 个，与 hk_genes_v2.txt 逐序一致）')
rec('C', f'JSON 内基因集 == 权威清单', summ.get('hk_genes_v2') == truth_genes,
    f"JSON {len(summ.get('hk_genes_v2',[]))} 个 / 清单 {len(truth_genes)} 个")
rows = re.findall(r"^\|\s*(\d+)\s*\|\s*([A-Z0-9]+)\s*\|\s*(\d+)\s*\|\s*([A-Z0-9]+)\s*\|\s*(\d+)\s*\|\s*([A-Z0-9]+)\s*\|$",
                  md_sum, re.M)
pairs = {}
for r in rows:
    pairs[int(r[0])] = r[1]; pairs[int(r[2])] = r[3]; pairs[int(r[4])] = r[5]
md_genes = [pairs[i] for i in range(1, 31)] if len(pairs) >= 30 else []
rec('C', '交付汇总 MD 基因表 == 权威清单', md_genes == truth_genes,
    f'MD 抽得 {len(md_genes)} 个')
rec('C', '结论差异 MD 含完整基因清单', truth_genes[0] in md_diff and truth_genes[-1] in md_diff)

# C2 富集率：HK 用显式计数，其余三组用区间（两者都是手稿的正确写法）
print('\n  C2. 逐表型富集率（GTEx ACAT-O）')
exp = {}
for r in summ['enrichment_acat_o']:
    exp[(r['Group'], r['Trait'])] = (r['N_FDR'], r['N_tested'], r['FDR_pct'])

print('    · HK 三组显式计数（手稿逐表型列出）')
for (g, t), (nf, nt, pct) in sorted(exp.items()):
    if g != 'Housekeeping v2':
        continue
    ok = (f'{nf}/{nt}' in en_text) and (f'{nf}/{nt}' in zh_text)
    rec('C', f'HK {t:<4} {nf}/{nt} = {pct}%  英+中均出现', ok)

print('    · 其余三组区间（手稿用区间表述，由 JSON 重算 min–max 校验）')
RANGE_STR = {'Candidate': '53.8–61.5', 'Non-Candidate': '46.2–69.2',
             'T2DM Control': '42.1–73.7'}
for g, expect in RANGE_STR.items():
    vals = [v[2] for (gg, t), v in exp.items() if gg == g]
    lo, hi = min(vals), max(vals)
    derived = f'{lo}–{hi}'
    ok_json = (derived == expect)
    ok_doc = (expect in en_text) and (expect in zh_text)
    rec('C', f'{g:<16} 区间 {expect}%  == JSON min–max 且英+中均出现',
        ok_json and ok_doc, f'JSON 逐表型 {sorted(vals)}')

# C3 合并富集率：由 JSON 重新累加，独立校验算术（防回归）
print('\n  C3. 三表型合并富集率（由 JSON 逐表型累加独立重算）')
pooled = {}
for (g, t), (nf, nt, pct) in exp.items():
    a, b = pooled.get(g, (0, 0))
    pooled[g] = (a + nf, b + nt)

EXPECT_POOLED = {'Housekeeping v2': (46, 84), 'Candidate': (46, 78),
                 'Non-Candidate': (43, 78), 'T2DM Control': (31, 57)}
for g in sorted(pooled):
    nf, nt = pooled[g]
    pct = round(100 * nf / nt, 1)
    in_docs = (str(pct) in en_text) and (str(pct) in zh_text)
    in_mds = (str(pct) in md_sum) and (str(pct) in md_diff)
    ok_counts = EXPECT_POOLED[g] == (nf, nt)
    # 手稿对三组只给区间、不给合并点估计（55.1/54.4 属附表细节），故只要求 MD 必含
    ok = ok_counts and in_mds and (in_docs or g in ('Non-Candidate', 'T2DM Control'))
    note = f'{nf}/{nt}={pct}% | 手稿{"含" if in_docs else "不含(区间表述)"} | 两份MD含'
    rec('C', f'{g:<16} 合并 {nf}/{nt} = {pct}%', ok, note)

rec('C', '合并分子/分母与交付值完全吻合（防回归）',
    all(pooled[g] == EXPECT_POOLED[g] for g in EXPECT_POOLED),
    ' | '.join(f'{g} {pooled[g][0]}/{pooled[g][1]}' for g in sorted(EXPECT_POOLED)))

# C4 Fisher
print('\n  C4. Fisher 精确检验（HK vs 各组，全部不显著）')
for p in ['0.64', '1.00']:
    ok = (p in en_text) and (p in zh_text)
    rec('C', f'Fisher P = {p} 英+中均出现', ok)

# C5 Figure 9
print('\n  C5. Figure 9 相关系数')
fig9 = summ['figure9']
for panel in ['a', 'b']:
    rho = fig9[panel]['rho']; n = fig9[panel]['n']
    s = f'{rho:.2f}'
    ok = (s in en_text) and (s in zh_text)
    rec('C', f'面板 ({panel}) ρ = {s}  n = {n}  英+中一致', ok,
        f'JSON ρ={rho:.6f}, n={n}')


# ============================================================ D. 无残留
print('\n' + '=' * 78)
print('D. 无残留检查（v1 污染披露 / 前稿表述）')
print('=' * 78)
PATS = ['earlier version', 'contaminated', 'superseded', 'supersedes', 'present revision',
        'reported previously', 'revised Table', 'revised Figure', '(v2)', 'discarded',
        '早期版本', '曾被污染', '取代了该早期', '本修订', '先前报告', '修订后的', '被弃用',
        '75.0', '15/20', 'violating its own']
for label, text in [('英文版', en_text), ('中文版', zh_text)]:
    hits = [(i, k) for i, p in enumerate(
        (doc_en if label == '英文版' else doc_zh).paragraphs)
        for k in PATS if k in p.text]
    if hits:
        rec('D', f'{label} 无前稿/污染表述', False, f'{len(hits)} 处命中: {hits[:5]}')
    else:
        rec('D', f'{label} 无前稿/污染表述', True, f'{len(PATS)} 个模式词全部零命中')


# ============================================================ 汇总
print('\n' + '=' * 78)
print('校验汇总')
print('=' * 78)
from collections import Counter
c = Counter((s, ok) for s, _, ok, _ in RESULTS)
for sec in ['A', 'B', 'C', 'D']:
    p, f = c[(sec, True)], c[(sec, False)]
    print(f'  {sec}: PASS {p} / FAIL {f}')
total_p = sum(v for (s, ok), v in c.items() if ok)
total_f = sum(v for (s, ok), v in c.items() if not ok)
print(f'\n  合计: PASS {total_p} / FAIL {total_f}')
if total_f == 0:
    print('  ✅ 全部通过——所有修改已落盘且一致无残留。')
else:
    print('  ⚠ 存在未通过项，请查看上方 FAIL 明细。')
