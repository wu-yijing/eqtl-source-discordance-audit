# -*- coding: utf-8 -*-
"""P0 Step 6：从官方数据重绘 Figure 8（跨性状三基因集组织轴）"""
import os, csv, re, shutil
import numpy as np
from scipy import stats
from docx import Document
from docx.table import Table
from docx.oxml.ns import qn
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import paths_config as P  # 统一路径入口（2026-09-20）

AF = P.need(P.AF1, 'Additional file 1（从期刊补充材料下载后用 AF1_DOCX 指定）')
NEW = P.need(P.DATA_Z, '官方 MetaXcan Z 数据层')
OUT = P.OUT_MAIN
BK = os.path.join(OUT, '_backup_before_officialZ_redraw_20260917')
os.makedirs(BK, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.7,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
C_TEST, C_HK, C_SCZ = '#C0392B', '#1E8449', '#2471A3'

d = Document(AF)
# 2026-09-23：tabs[i] ≠ Table S(i+1)。实测本版 AF1 有 25 个表体，其中 **Note S1 自带一张表**，
#   使旧的硬编码索引整体错位 1（grid(1) 取到 Table S1、grid(5) 取到 Table S4b），脚本直接
#   IndexError：grid(1) 是 7 列的基因清单，没有 Z_Whole_Blood 列 → 无一行动过 len(r)>=10 的过滤 →
#   is_cand 变成空的 float64 数组 → 报 "arrays used as indices must be of integer (or boolean) type"。
#   改为按 body 顺序的题注解析表号（与 AF1 自身编号一致，加表/挪表都不会再错位）。
_by_label, _cur = {}, None
for _ch in d.element.body.iterchildren():
    if _ch.tag == qn('w:p'):
        _t = ''.join(_x.text or '' for _x in _ch.iter(qn('w:t'))).strip()
        _m = re.match(r'^\**\s*Table\s+(S\d+[a-c]?)\b', _t)
        if _m:
            _cur = _m.group(1)
    elif _ch.tag == qn('w:tbl'):
        if _cur and _cur not in _by_label:
            _by_label[_cur] = Table(_ch, d)
print('AF1 tables resolved:', sorted(_by_label))
grid = lambda s: [[c.text.strip() for c in r.cells] for r in _by_label[s].rows]
f = lambda x: (float(x) if x not in (None, '', 'NA') else np.nan)

def rho_ci(x, y):
    m = ~(np.isnan(x) | np.isnan(y))
    a, b = x[m], y[m]
    r = stats.spearmanr(a, b)[0]
    n = len(a)
    z = np.arctanh(r); se = 1 / np.sqrt(n - 3)
    return r, n, np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)

def rdc(x, y):
    m = ~(np.isnan(x) | np.isnan(y))
    return 100 * ((x[m] > 0) == (y[m] > 0)).mean()

# --- 测试集：tissue-only（GTEx WB vs NT），全部有 Z 的基因 ---
S2 = grid('S2'); i2 = {v: k for k, v in enumerate(S2[0])}
G = {r[0]: r[1] for r in grid('S1')[1:] if r[0] and 'andidate' in r[1] or (r[0] and r[1])}
zwb = np.array([f(r[i2['Z_Whole_Blood']]) for r in S2[1:] if len(r) >= 10 and r[0]])
znt = np.array([f(r[i2['Z_Nerve_Tibial']]) for r in S2[1:] if len(r) >= 10 and r[0]])
genes2 = [r[0] for r in S2[1:] if len(r) >= 10 and r[0]]
is_cand = np.array([('andidate' in str(G.get(g, '')) and 'Non' not in str(G.get(g, ''))) or
                    ('Non-Candidate' in str(G.get(g, ''))) or ('T2DM' in str(G.get(g, ''))) for g in genes2])
r_t, n_t, lo_t, hi_t = rho_ci(zwb[is_cand], znt[is_cand])
print('testbed tissue-only: rho=%.3f n=%d CI %.3f-%.3f  consistency=%.1f%%'
      % (r_t, n_t, lo_t, hi_t, rdc(zwb[is_cand], znt[is_cand])))

# --- housekeeping ---
S5 = grid('S5'); h5 = S5[0]
zi = [i for i, h in enumerate(h5) if 'Z, Nerve_Tibial' in h]
zj = [i for i, h in enumerate(h5) if 'Z, Whole_Blood' in h]
print('S5 cols NT:', zi, 'WB:', zj)
hk_nt, hk_wb = [], []
for r in S5[1:]:
    if len(r) <= max(zi + zj) or not r[0]:
        continue
    for a, b in zip(zi, zj):
        hk_nt.append(f(r[a])); hk_wb.append(f(r[b]))
hk_nt = np.array(hk_nt); hk_wb = np.array(hk_wb)
r_hk, n_hk, lo_hk, hi_hk = rho_ci(hk_wb, hk_nt)
print('housekeeping tissue-only: rho=%.3f n=%d CI %.3f-%.3f' % (r_hk, n_hk, lo_hk, hi_hk))

# --- SCZ genome-wide ---
rows = list(csv.DictReader(open(os.path.join(NEW, 'scz_z_4arm_official.csv'), encoding='utf-8-sig')))
g = lambda k: np.array([float(r[k]) if r[k] not in ('', 'NA') else np.nan for r in rows])
ze, zwb_s, znt_s = g('eqZ'), g('wbZ'), g('ntZ')
r_s, n_s, lo_s, hi_s = rho_ci(zwb_s, znt_s)
print('SCZ tissue-only: rho=%.3f n=%d CI %.3f-%.3f  consistency=%.1f%%'
      % (r_s, n_s, lo_s, hi_s, rdc(zwb_s, znt_s)))

panels = [('HOTAIR testbed\n(FinnGen R13, official Z)', r_t, n_t, lo_t, hi_t, C_TEST),
          ('Housekeeping control\n(FinnGen R13, official Z)', r_hk, n_hk, lo_hk, hi_hk, C_HK),
          ('PGC3 schizophrenia\n(genome-wide, official Z)', r_s, n_s, lo_s, hi_s, C_SCZ)]
fig, ax = plt.subplots(figsize=(5.2, 3.2))
xs = np.arange(3)
# 2026-09-23：首/末点的 ρ/n 标注改为向内对齐。原为三点一律 ha='center'，首点文字以 x=0 居中，
#   左半压到 y 轴脊线上（ρ 字形与脊线重叠，并与 0.6 刻度标签相挤）。ax.text 默认 clip_on=False，
#   所以是"压住"而非"被裁"。两侧点改为贴内侧对齐，且不改动 xlim —— 数据点与误差线位置逐像素不变。
_ALIGN = {0: ('left', 0.04), 1: ('center', 0.0), 2: ('right', -0.04)}
for i, (lab, r, n, lo, hi, c) in enumerate(panels):
    ax.errorbar(i, r, yerr=[[r - lo], [hi - r]], fmt='o', ms=6.5, color=c, ecolor=c, capsize=4, lw=1.2)
    _ha, _dx = _ALIGN[i]
    ax.text(i + _dx, hi + 0.03, '$\\rho$ = %+.2f\nn = %d' % (r, n), ha=_ha, fontsize=6.2)
ax.set_xticks(xs); ax.set_xticklabels([p[0] for p in panels], fontsize=6.4)
ax.set_ylabel("Spearman $\\rho$ (GTEx Whole_Blood vs Nerve_Tibial)")
_hi = max(hi for _, _, _, _, hi, _ in panels)
ax.set_ylim(0.30, max(0.80, _hi + 0.11))
ax.axhline(0, color='#AAAAAA', lw=0.5)
# 2026-09-23：移除图内标题。BMC 要求标题只出现在图注里；该标题此前是**手工从输出图上抹掉**的，
#   脚本一直没同步 —— 直接重跑会把标题带回来。此处移除，使脚本产出与定稿图一致。


def save(fig, name):
    """2026-09-23：补上 RGBA→RGB 转换。01_redraw_Fig5_Fig7.py 在 2026-09-20 已加，本脚本当时漏了；
    不补的话重跑会用 RGBA 覆盖定稿图里的 RGB 版本（图集其余图均为 RGB）。"""
    for ext in ['png', 'pdf']:
        src = os.path.join(OUT, f'{name}.{ext}')
        if os.path.exists(src) and not os.path.exists(os.path.join(BK, f'{name}.{ext}')):
            shutil.copy2(src, os.path.join(BK, f'{name}.{ext}'))
        fig.savefig(src, dpi=600)
    from PIL import Image as _Image
    _p = os.path.join(OUT, name + '.png')
    _im = _Image.open(_p)
    if _im.mode != 'RGB':
        _im.convert('RGB').save(_p, dpi=(600, 600))
    plt.close(fig)


fig.tight_layout(); save(fig, 'Fig8')
print('Fig8 redrawn')
print('backups:', sorted(os.listdir(BK)))
