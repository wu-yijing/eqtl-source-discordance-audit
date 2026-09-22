# -*- coding: utf-8 -*-
"""Fig. S6 重出（2026-09-20）。

修复
----
原脚本（`2026-09-17-20-59-12/_review/apply_edits_2.py`）把 BH 检测边界 |Z|（实测 2.97-3.44）
与单基因 power（0-1）画在同一根 `axa` 上，而 `axa.set_ylim(0, 1.18)` → 边界点全部落到坐标轴之外，
面板 (a) 只剩坐标轴/刻度/图例（图例首项 "BH boundary |Z| at q = 0.05" 的实测值 2.97-3.69 与
0-1.1 的量纲自相矛盾，即此因）。面板 (b) 图例 `loc='upper left'` 落在绘图区内，压住功率曲线。

修法
----
* 面板 (a)：双 y 轴 —— 左轴 = 检测边界 |Z|（对数 stratum size），右轴 = 单基因 power（k = 1）。
* 面板 (b)：图例移到绘图区之外（轴下方），不再压线。
* 内置回归断言：写图的每一个 n 的检测边界必须与 SI Table S19 的 "Detection boundary |Z|" 一致。
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import paths_config as P  # 统一路径入口（2026-09-20）

WD = P.RES
OUT = P.OUT_SUPP

# SI Table S19 的检测边界（逐行抄自 Additional file 1，2026-09-20 核验）
S19 = {28: 3.12, 19: 3.01, 29: 3.13, 87: 3.44, 84: 3.43,
       27: 3.11, 17: 2.97, 81: 3.42}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.7,
                     'axes.spines.top': False, 'pdf.fonttype': 42, 'ps.fonttype': 42})
C1, C2, C3 = '#C0392B', '#2471A3', '#7D3C98'


def main():
    PC = json.load(open(os.path.join(WD, 'm15_positive_control.json'), encoding='utf-8'))

    # ---- 回归断言：边界值与 SI Table S19 逐点一致 ----
    bad = []
    for n, v in S19.items():
        got = PC['PC1a_BH_boundary'].get(str(n), {}).get('min_absZ')
        if got is None or abs(float(got) - v) > 5e-3:
            bad.append((n, got, v))
    if bad:
        print('!! 检测边界与 Table S19 不符:', bad)
        return 1
    print('OK 检测边界与 Table S19 的 %d 个层逐点一致' % len(S19))

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(6.65, 3.0))

    # ================= 面板 (a)：双 y 轴 =================
    nb = sorted(int(k) for k in PC['PC1a_BH_boundary'] if int(k) <= 300)
    bvals = [float(PC['PC1a_BH_boundary'][str(n)]['min_absZ']) for n in nb]
    axa.semilogx(nb, bvals, 'o-', color=C1, ms=3.4, lw=1.1,
                 label=r'BH boundary $|Z|$ at $q=0.05$')
    axa.axhline(3.0, color='#888888', lw=0.7, ls=':')
    axa.set_xlim(14, 300)
    lo, hi = min(bvals), max(bvals)
    axa.set_ylim(lo - 0.35, hi + 0.35)
    axa.set_ylabel(r'BH detection boundary $|Z|$', fontsize=7.0)
    axa.set_xlabel(r'stratum size $n$', fontsize=7.0)
    axa.tick_params(axis='both', labelsize=6.4)
    # 2026-09-20：log 轴的 10^n 刻度会把指数渲染成 0.7x 字号（6.4 -> 4.48 pt，低于 BMC 6 pt 下限），
    # 故改用显式刻度标签，彻底消除上标。
    axa.set_xticks([20, 50, 100, 200])
    axa.set_xticklabels(['20', '50', '100', '200'])

    axr = axa.twinx()
    for n, c, mk in ((17, C3, 's'), (27, C2, '^'), (84, '#1E8449', 'D')):
        lab = [k for k in PC['PC2a_single_gene_power'] if '(n=%d)' % n in k]
        if not lab:
            continue
        d = PC['PC2a_single_gene_power'][lab[0]]
        lams = sorted(set(float(k.split('lam=')[1]) for k in d['curves'] if k.startswith('k=1|')))
        ys = [d['curves']['k=1|lam=%s' % l] for l in lams]
        axr.plot(lams, ys, mk + '--', ms=2.8, lw=1.0, color=c, alpha=0.9,
                 label=r'power, $k=1$, $n=%d$' % n)
    axr.axhline(0.8, color='#555555', lw=0.7, ls=':')
    axr.set_ylim(0, 1.18)
    axr.set_ylabel('single-gene power', fontsize=7.0)
    axr.tick_params(axis='y', labelsize=6.4)
    axr.spines['top'].set_visible(False)

    h1, l1 = axa.get_legend_handles_labels()
    h2, l2 = axr.get_legend_handles_labels()
    axa.legend(h1 + h2, l1 + l2, frameon=False, fontsize=6.2, loc='lower right',
               handlelength=1.8, borderpad=0.2)
    axa.set_title('(a) Endpoint detection boundary and single-gene power', fontsize=7.6, pad=4)

    # ================= 面板 (b)：图例移出绘图区 =================
    pb = PC['PC2b_group_diff_power']
    scen = [('GTEx: HK 0/87 vs cand 2/84', pb['GTEx (HK 87 @0.0% vs cand 84 @2.4%)'], C1, 'o'),
            ('eQTLGen: HK 2/81 vs cand 5/81', pb['eQTLGen (HK 81 @2.5% vs cand 81 @6.2%)'], C2, 's'),
            ('GTEx nominal 5.7% vs 9.5%', pb['GTEx (HK 87 @5.7% vs cand 84 @9.5%)'], C1, '^'),
            ('eQTLGen nominal 8.6% vs 9.9%', pb['eQTLGen (HK 81 @8.6% vs cand 81 @9.9%)'], C2, 'v')]
    for lab, d, c, mk in scen:
        xs = sorted(int(k) for k in d['power_curve'])
        axb.plot(xs, [d['power_curve'][str(x)] for x in xs], marker=mk, ms=3.2, lw=1.05,
                 color=c, label=lab)
    axb.axhline(0.8, color='#555555', lw=0.7, ls=':')
    axb.text(0.4, 0.83, '80% power', fontsize=6.0, color='#555555')
    axb.axvline(2.7, color=C3, lw=0.8, ls='--')
    axb.text(2.95, 0.62, 'observed\n+2.7 pp', fontsize=6.0, color=C3, ha='left', va='center')
    axb.set_xlabel('true candidate-over-housekeeping difference (pp)', fontsize=7.0)
    axb.set_ylabel('power (two-sided exact Fisher)', fontsize=7.0)
    axb.set_ylim(0, 1.02)
    axb.set_xlim(0, 20)
    axb.tick_params(axis='both', labelsize=6.4)
    axb.legend(frameon=False, fontsize=6.2, loc='upper center',
               bbox_to_anchor=(0.5, -0.30), ncol=2, handlelength=1.8, columnspacing=1.0)
    axb.set_title('(b) Between-group power of the enrichment endpoint', fontsize=7.6, pad=4)
    axb.spines['top'].set_visible(False)

    fig.tight_layout(rect=(0, 0.06, 1, 1))
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, 'FigS6.%s' % ext), dpi=600, facecolor='white')
    plt.close(fig)

    # PNG 转 RGB
    p = os.path.join(OUT, 'FigS6.png')
    im = Image.open(p)
    if im.mode != 'RGB':
        im.convert('RGB').save(p, dpi=(600, 600))
    im = Image.open(p)
    print('OK FigS6 written: %s mode=%s size=%s dpi=%s' % (p, im.mode, im.size, im.info.get('dpi')))
    print('   pdf bytes =', os.path.getsize(os.path.join(OUT, 'FigS6.pdf')))
    return 0


if __name__ == '__main__':
    sys.exit(main())
