# -*- coding: utf-8 -*-
"""P0 Step 7：从官方数据重算三个分解臂，判定 Fig.4 / Table 3b / 正文三套值哪个正确"""
import os, csv
import numpy as np
from scipy import stats
import paths_config as P  # 统一路径入口（2026-09-20）

NEW = P.need(P.DATA_Z, '官方 MetaXcan Z 数据层')
def rd(n):
    return list(csv.DictReader(open(os.path.join(NEW, n), encoding='utf-8-sig')))
f = lambda x: (float(x) if x not in (None, '', 'NA') else np.nan)

GT = {(r['Gene'], r['Trait']): r for r in rd('gtex_official_Z.csv')}
EQ = {(r['Gene'], r['Trait']): r for r in rd('eqtlgen_official_Z.csv')}
GRP = {r['Gene']: r['Group'] for r in rd('gene_groups_TableS1_official.csv')}
keys = [k for k in EQ]
print('eQTLGen 测试集对数:', len(keys))

def arm(name, fa, fb, universe):
    """fa/fb: 取 Z 的函数"""
    xs, ys, tr = [], [], []
    for k in universe:
        a, b = fa(k), fb(k)
        if np.isnan(a) or np.isnan(b):
            continue
        xs.append(a); ys.append(b); tr.append(k[1])
    xs, ys = np.array(xs), np.array(ys)
    if len(xs) < 5:
        return None
    r = stats.spearmanr(xs, ys)[0]
    return dict(arm=name, n=len(xs), rho=round(float(r), 4),
                consistency=round(100 * float(((xs > 0) == (ys > 0)).mean()), 1),
                k=int(((xs > 0) == (ys > 0)).sum()))

Zg_wb = lambda k: f(GT[k]['Z_Whole_Blood']) if k in GT else np.nan
Zg_nt = lambda k: f(GT[k]['Z_Nerve_Tibial']) if k in GT else np.nan
Zg_mt = lambda k: f(GT[k]['Z_multi_tissue']) if k in GT else np.nan
Ze = lambda k: f(EQ[k]['Z_eQTLGen']) if k in EQ else np.nan

# 宇宙定义
allk = sorted(set(GT) | set(EQ))
pan = [k for k in allk if not np.isnan(Ze(k)) and not np.isnan(Zg_wb(k))]
tis = [k for k in allk if not np.isnan(Zg_wb(k)) and not np.isnan(Zg_nt(k))]
dua = [k for k in allk if not np.isnan(Ze(k)) and not np.isnan(Zg_mt(k))]
common = [k for k in allk if not any(np.isnan(v) for v in
          (Ze(k), Zg_wb(k), Zg_nt(k), Zg_mt(k)))]
print('panel-only 宇宙 %d 对 / %d 基因' % (len(pan), len(set(x[0] for x in pan))))
print('tissue-only 宇宙 %d 对 / %d 基因' % (len(tis), len(set(x[0] for x in tis))))
print('dual 宇宙 %d 对 / %d 基因' % (len(dua), len(set(x[0] for x in dua))))
print('四臂共同宇宙 %d 对 / %d 基因' % (len(common), len(set(x[0] for x in common))))
gpan = sorted(set(x[0] for x in pan)); gtis = sorted(set(x[0] for x in tis)); gdua = sorted(set(x[0] for x in dua))
gcom = sorted(set(x[0] for x in common))
print('面板基因数 %d / %d / %d ; 共同 %d' % (len(gpan), len(gtis), len(gdua), len(gcom)))

print('\n=== (A) 各臂在自身完整宇宙（full-pair）===')
for n, fa, fb, u in [('panel-only', Ze, Zg_wb, pan), ('tissue-only', Zg_wb, Zg_nt, tis),
                     ('dual', Ze, Zg_mt, dua)]:
    print('  ', arm(n, fa, fb, u))
print('=== (B) 四臂共同宇宙 ===')
for n, fa, fb in [('panel-only', Ze, Zg_wb), ('tissue-only', Zg_wb, Zg_nt), ('dual', Ze, Zg_mt)]:
    print('  ', arm(n, fa, fb, common))
print('\n=== 稿件中的三套值 ===')
print('  Table 3b full-pair : panel 0.490 (71.1%) / tissue 0.414 (65.9%) / dual 0.442 (67.2%)')
print('  正文 §3.3 45 基因共同宇宙: panel 0.463 (70.4%) / tissue 0.410 (65.2%) / dual 0.443 (67.4%)')
print('  Fig.4 脚本硬编码     : panel 0.355 (67.3%) / tissue 0.422 (67.3%) / dual 0.287 (59.2%)')
