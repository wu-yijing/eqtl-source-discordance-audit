# -*- coding: utf-8 -*-
"""P0 Step 5：从官方 Z 重绘 Figure 3（头条散点图）与 Figure 6（富集三面板）"""
import os, csv, shutil
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import paths_config as P  # 统一路径入口（2026-09-20）

NEW = P.need(P.DATA_Z, '官方 MetaXcan Z 数据层')
OUT = P.OUT_MAIN
BK = os.path.join(OUT, '_backup_before_officialZ_redraw_20260917')
os.makedirs(BK, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.7,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
C_GTEX, C_EQTL = '#C0392B', '#2471A3'
CGRP = {'Candidate': '#C0392B', 'Non-Candidate': '#2471A3', 'T2DM control': '#1E8449'}

def rd(name):
    with open(os.path.join(NEW, name), encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))
def save(fig, name):
    # 2026-09-20 拆分：Fig6 已由 08_redraw_Fig6_labels_20260920.py 接管，
    # 本脚本的 Fig6 段是旧版（ylim=24、图例在内），重跑会覆盖正确版，故默认跳过。
    if name == 'Fig6' and not os.environ.get('FIG6_FROM_02'):
        print('SKIP Fig6 —— 权威脚本是 08_redraw_Fig6_labels_20260920.py')
        plt.close(fig)
        return
    for ext in ['png', 'pdf']:
        src = os.path.join(OUT, f'{name}.{ext}')
        if os.path.exists(src) and not os.path.exists(os.path.join(BK, f'{name}.{ext}')):
            shutil.copy2(src, os.path.join(BK, f'{name}.{ext}'))
        fig.savefig(src, dpi=600)
    # 2026-09-20：matplotlib 默认写 RGBA，而图集其余图为 RGB（评审 m-14），统一转 RGB。
    from PIL import Image as _Image
    _p = os.path.join(OUT, name + '.png')
    _im = _Image.open(_p)
    if _im.mode != 'RGB':
        _im.convert('RGB').save(_p, dpi=(600, 600))
    plt.close(fig)
f = lambda x: (float(x) if x not in (None, '', 'NA') else np.nan)

# ---------------- Figure 3 ----------------
S12 = rd('primary_arm_96pairs_official.csv')
zg = np.array([f(r['Z_GTEx']) for r in S12]); ze = np.array([f(r['Z_eQTLGen']) for r in S12])
tr = np.array([r['Trait'] for r in S12])
same = (zg > 0) == (ze > 0)
rho, pnaive = stats.spearmanr(zg, ze)
print('Fig3: n=%d consistency=%.1f%% rho=%.4f naiveP=%.2g max|Ze|=%.2f max|Zg|=%.2f'
      % (len(S12), 100 * same.mean(), rho, pnaive, np.abs(ze).max(), np.abs(zg).max()))
fig, (axa, axb) = plt.subplots(1, 2, figsize=(6.65, 3.1))
cols = [C_GTEX if s else '#7F8C8D' for s in same]
axa.scatter(zg, ze, s=14, c=cols, edgecolor='none', alpha=0.85, zorder=3)
lim = max(np.abs(zg).max(), np.abs(ze).max()) * 1.18
axa.plot([-lim, lim], [-lim, lim], ls='--', lw=0.8, color='#555555', zorder=2)
axa.axhline(0, color='#AAAAAA', lw=0.5); axa.axvline(0, color='#AAAAAA', lw=0.5)
axa.set_xlim(-lim, lim); axa.set_ylim(-lim, lim)
axa.set_xlabel('GTEx v8 Whole_Blood Z'); axa.set_ylabel('eQTLGen whole-blood Z')
axa.text(0.03, 0.96, 'Spearman $\\rho$ = %.2f\nDirection consistency = %.1f%% (%d/%d)\nnaive P = %.1g (anticonservative)'
         % (rho, 100 * same.mean(), int(same.sum()), len(same), pnaive),
         transform=axa.transAxes, va='top', ha='left', fontsize=6.2)
axa.set_title('(a) 96 primary-arm gene\u2013phenotype pairs (official Z)', fontsize=7.6, pad=4)
order = ['DR', 'DN', 'DPN']
data = [zg[tr == t] for t in order]
bp = axb.boxplot(data, positions=range(3), widths=0.55, patch_artist=True, showfliers=True,
                 flierprops=dict(marker='o', ms=2.4, alpha=0.6, markerfacecolor=C_GTEX, markeredgecolor='none'))
for pc in bp['boxes']:
    pc.set_facecolor(C_GTEX); pc.set_alpha(0.35)
for med in bp['medians']:
    med.set_color('#1A1A1A'); med.set_linewidth(1.1)
axb.set_ylim(top=axb.get_ylim()[1] * 1.20)   # 顶部留白：避免 consistency 标注压到面板标题（2026-09-20）
for i, t in enumerate(order):
    axb.text(i, 0.965, 'consistency %.1f%%' % (100 * same[tr == t].mean()), ha='center', va='top',
             fontsize=6.3, transform=axb.get_xaxis_transform())
axb.axhline(0, color='#AAAAAA', lw=0.5)
axb.set_xticks(range(3)); axb.set_xticklabels(order)
axb.set_ylabel('GTEx v8 Whole_Blood Z')
axb.set_title('(b) GTEx Z distribution by phenotype', fontsize=7.6, pad=4)
fig.tight_layout(); save(fig, 'Fig3')
print('Fig3 redrawn')

# ---------------- Figure 6 ----------------
GT = rd('gtex_official_Z.csv'); EQ = rd('eqtlgen_official_Z.csv')
GRP = {r['Gene']: r['Group'] for r in rd('gene_groups_TableS1_official.csv')}
def norm(g):
    return ('Candidate' if ('Candidate' in g and 'Non' not in g) else
            'Non-Candidate' if 'Non-Candidate' in g else
            'T2DM control' if 'T2DM' in g else g)
GTX, EQX = {}, {}
for r in rd('gtex_official_Z.csv'):
    p = f(r['P_ACAT_O'])
    g = norm(GRP.get(r['Gene'], ''))
    if not np.isnan(p):
        GTX.setdefault((g, r['Trait']), []).append(p)
for r in rd('eqtlgen_official_Z.csv'):
    p = f(r['P'])
    if not np.isnan(p):
        EQX.setdefault((norm(r['Group']), r['Trait']), []).append(p)
GRPS = ['Candidate', 'Non-Candidate', 'T2DM control', 'Housekeeping']
PH = ['DR', 'DN', 'DPN']
def cp(k, n, a=0.05):
    lo = 0.0 if k == 0 else stats.beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return 100 * lo, 100 * hi
fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0), sharey=True)
x = np.arange(len(GRPS)); w = 0.36
for pi, ph in enumerate(PH):
    ax = axes[pi]
    for si, (src, col, off) in enumerate([('GTEx v8 multi-tissue ACAT-O', C_GTEX, -w / 2),
                                          ('eQTLGen whole blood', C_EQTL, w / 2)]):
        vals, los, his, ns = [], [], [], []
        for g in GRPS:
            ps = np.array((GTX if src.startswith("GTEx") else EQX).get((g, ph), []))
            if len(ps) == 0:
                vals.append(0); los.append(0); his.append(0); ns.append((0, 0)); continue
            k = int((ps < 0.05).sum()); n = len(ps)
            r = 100 * k / n; lo, hi = cp(k, n)
            vals.append(r); los.append(r - lo); his.append(hi - r); ns.append((k, n))
        ax.bar(x + off, vals, w, color=col, edgecolor='black', lw=0.5,
               label=src if pi == 0 else None)
        ax.errorbar(x + off, vals, yerr=[los, his], fmt='none', ecolor='#222222', elinewidth=0.8, capsize=3)
        for j, v in enumerate(vals):
            ax.text(x[j] + off, v + his[j] + 0.45, '%d/%d' % ns[j], ha='center', fontsize=5.2)
    ax.set_xticks(x); ax.set_xticklabels(['Cand', 'Non-cand', 'T2DM', 'HK'], fontsize=6.4)
    ax.set_ylim(0, 24); ax.set_title('(%s) %s' % ('abc'[pi], ph), fontsize=7.6, pad=4)
    if pi == 0:
        ax.set_ylabel('Genes with nominal p < 0.05 (%)')
        ax.legend(frameon=False, fontsize=6.0, loc='upper left')
fig.tight_layout(); save(fig, 'Fig6')
print('Fig6 redrawn')
for g in GRPS:
    row = []
    for ph in PH:
        for src in ['GTEx', 'eQTLGen']:
            pass
    for src, DD in [('GTEx', GTX), ('eQTLGen', EQX)]:
        print('   %-14s %-12s' % (src, g), {ph: '%d/%d' % (int((np.array(DD.get((g, ph), [0])) < 0.05).sum()), len(DD.get((g, ph), []))) for ph in PH})
print('\nbackups ->', BK, sorted(os.listdir(BK)))
