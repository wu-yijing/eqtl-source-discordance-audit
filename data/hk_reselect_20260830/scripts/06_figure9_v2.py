#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 9 v2：三面板。a/b 为看家基因对照内 GTEx 对比；c 为组织语境轴的跨基因集比较。"""
import io, sys, os, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
ROOT = r'E:\workbuddy\hk_reselect_20260830'


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    rx, ry = pd.Series(x).rank().values, pd.Series(y).rank().values
    mx, my = rx.mean(), ry.mean()
    den = math.sqrt(((rx - mx) ** 2).sum() * ((ry - my) ** 2).sum())
    rho = ((rx - mx) * (ry - my)).sum() / den if den > 0 else float('nan')
    z = 0.5 * math.log((1 + rho) / (1 - rho))
    se = 1 / math.sqrt(n - 3)
    return rho, math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se), n


arms = pd.read_csv(os.path.join(ROOT, 'arms_all_groups.csv'))
h = arms[arms.grp == 'HK_v2']

fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8))

specs = [('z_ACAT_O', 'z_Nerve_Tibial', 'GTEx multi-tissue (ACAT-O) Z', 'GTEx Nerve_Tibial Z', 'a'),
         ('z_Whole_Blood', 'z_Nerve_Tibial', 'GTEx Whole_Blood Z', 'GTEx Nerve_Tibial Z', 'b')]
stats = {}
for ax, (xc, yc, xl, yl, lab) in zip(axes[:2], specs):
    s = h.dropna(subset=[xc, yc])
    rho, lo, hi, n = spearman(s[xc].values, s[yc].values)
    stats[lab] = (rho, lo, hi, n)
    ax.scatter(s[xc], s[yc], s=38, facecolors='none', edgecolors='#1f77b4', linewidths=1.3)
    lim = max(np.abs(s[xc]).max(), np.abs(s[yc]).max()) * 1.15
    ax.plot([-lim, lim], [-lim, lim], ls='--', lw=1, color='#888888')
    ax.axhline(0, lw=0.6, color='#cccccc')
    ax.axvline(0, lw=0.6, color='#cccccc')
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel(xl, fontsize=10); ax.set_ylabel(yl, fontsize=10)
    ax.tick_params(labelsize=9)
    ax.text(0.03, 0.96, f'Spearman rho = {rho:+.2f}\n95% CI {lo:+.2f} to {hi:+.2f}\nn = {n} pairs',
            transform=ax.transAxes, ha='left', va='top', fontsize=9.5)
    ax.text(0.5, -0.15, f'({lab}) Housekeeping true-negative control (v2)',
            transform=ax.transAxes, ha='center', va='top', fontsize=10, fontweight='bold')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

# ---- 面板 c：组织语境轴跨基因集比较 ----
ax = axes[2]
labels = ['PGC SCZ3\n(pooled)', 'HOTAIR\ntestbed', 'Housekeeping\nv2']
vals = [0.51, 0.45, stats['b'][0]]
ns = [2511, 150, stats['b'][3]]
err = [[0, 0, vals[2] - stats['b'][1]], [0, 0, stats['b'][2] - vals[2]]]
colors = ['#9ecae1', '#a1d99b', '#fdae6b']
bars = ax.bar(range(3), vals, color=colors, edgecolor='#555555', linewidth=0.8, width=0.6)
ax.errorbar(range(3), vals, yerr=err, fmt='none', ecolor='#333333', capsize=5, lw=1.2)
for i, (v, n) in enumerate(zip(vals, ns)):
    ax.text(i, v + 0.035, f'{v:+.2f}\nn = {n}', ha='center', va='bottom', fontsize=9.5)
ax.set_xticks(range(3))
ax.set_xticklabels(labels, fontsize=9.5)
ax.set_ylabel('Tissue-context axis: Spearman rho\n(GTEx Whole_Blood vs GTEx Nerve_Tibial)', fontsize=9.5)
ax.set_ylim(0, 0.95)
ax.tick_params(labelsize=9)
ax.text(0.5, -0.15, '(c) Tissue-context axis across gene sets',
        transform=ax.transAxes, ha='center', va='top', fontsize=10, fontweight='bold')
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)

plt.tight_layout()
plt.subplots_adjust(bottom=0.24, wspace=0.30)
plt.savefig(os.path.join(ROOT, 'Figure9_hk_v2.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(ROOT, 'Figure9_hk_v2.pdf'), bbox_inches='tight')

print('Figure 9 v2 统计:')
print(f'  (a) multi-tissue vs Nerve_Tibial : rho = {stats["a"][0]:+.3f} '
      f'95%CI [{stats["a"][1]:+.2f}, {stats["a"][2]:+.2f}]  n = {stats["a"][3]}')
print(f'  (b) Whole_Blood vs Nerve_Tibial  : rho = {stats["b"][0]:+.3f} '
      f'95%CI [{stats["b"][1]:+.2f}, {stats["b"][2]:+.2f}]  n = {stats["b"][3]}')
print(f'  (c) 组织语境轴: PGC 0.51 (n=2511) / HOTAIR 0.45 (n=150) / HK {stats["b"][0]:+.2f} (n={stats["b"][3]})')
