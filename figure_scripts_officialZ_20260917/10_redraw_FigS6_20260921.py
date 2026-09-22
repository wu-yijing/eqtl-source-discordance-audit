# -*- coding: utf-8 -*-
"""Fig. S6 重出（2026-09-21）—— 承接 `10_redraw_FigS6_20260920.py`。

本轮修 4 项（起因：图注合规审计，2026-09-21）
--------------------------------------------
1. **x 轴补标 log**：xlabel 由 `stratum size n` 改为 `stratum size n (log scale)`。
   （`semilogx` 一直是对数轴，但轴标签没写，读者会把 20/50/100/200 的间距误读成线性。）
2. **两条参考线补入图内图例**（BMC：「Figure keys should be incorporated into the graphic,
   not into the legend of the figure.」）：
   - `axa.axhline(3.0)` → 图例项 `|Z| = 3.0`（左轴）
   - `axr.axhline(0.8)` → 图例项 `power = 0.8`（右轴）
   两条线原本在整个图里没有任何文字说明（面板 (b) 的同类线都有标注）。
3. **数学符号改正体**：mathtext 默认把变量排成斜体（DejaVuSans-Oblique），
   而全稿 house style 是正体 —— 实测主稿 118 处 + SI 83 处含统计符号的 run，
   斜体数 = **0**（主稿 8 条图注只斜体基因名）。故用 `\\mathrm{}` 把 Z/q/k/n 改为正体，
   使图件与图注、与其余 13 条图注一致。
4. 保留 2026-09-20 版的全部修法（双 y 轴、图例移出绘图区、与 Table S19 的回归断言）。

注意：`S19` 断言表只覆盖 Table S19 中**已列表**的 8 个层，而 `PC1a_BH_boundary` 里
n ≤ 300 的点全都会被画出（含 n = 222 = GTEx multi-tissue ACAT-O 主臂的 222 个
gene–phenotype pairs）。Table S19 已于 2026-09-21 补入 n = 222 一行，本脚本不改表。
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
       27: 3.11, 17: 2.97, 81: 3.42, 222: 3.69}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.7,
                     'axes.spines.top': False, 'pdf.fonttype': 42, 'ps.fonttype': 42})
C1, C2, C3 = '#C0392B', '#2471A3', '#7D3C98'

L_BH = r'BH boundary $|\mathrm{Z}|$ at $\mathrm{q}=0.05$'
L_Z3 = r'$|\mathrm{Z}| = 3.0$'
L_P8 = r'power $= 0.8$'


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
    axa.semilogx(nb, bvals, 'o-', color=C1, ms=3.4, lw=1.1, label=L_BH)
    axa.axhline(3.0, color='#888888', lw=0.7, ls=':', label=L_Z3)
    axa.set_xlim(14, 300)
    lo, hi = min(bvals), max(bvals)
    axa.set_ylim(lo - 0.35, hi + 0.35)
    axa.set_ylabel(r'BH detection boundary $|\mathrm{Z}|$', fontsize=7.0)
    axa.set_xlabel(r'stratum size $\mathrm{n}$ (log scale)', fontsize=7.0)
    axa.tick_params(axis='both', labelsize=6.4)
    # 2026-09-20：log 轴的 10^n 刻度会把指数渲染成 0.7x 字号（6.4 -> 4.48 pt，低于 BMC 6 pt 下限），
    # 故改用显式刻度标签，彻底消除上标。
    axa.set_xticks([20, 50, 100, 200])
    axa.set_xticklabels(['20', '50', '100', '200'])

    axr = axa.twinx()
    power_labs = []
    for n, c, mk in ((17, C3, 's'), (27, C2, '^'), (84, '#1E8449', 'D')):
        lab = [k for k in PC['PC2a_single_gene_power'] if '(n=%d)' % n in k]
        if not lab:
            continue
        d = PC['PC2a_single_gene_power'][lab[0]]
        lams = sorted(set(float(k.split('lam=')[1]) for k in d['curves'] if k.startswith('k=1|')))
        ys = [d['curves']['k=1|lam=%s' % l] for l in lams]
        pl = r'power, $\mathrm{k}=1$, $\mathrm{n}=%d$' % n
        axr.plot(lams, ys, mk + '--', ms=2.8, lw=1.0, color=c, alpha=0.9, label=pl)
        power_labs.append(pl)
    axr.axhline(0.8, color='#555555', lw=0.7, ls=':', label=L_P8)
    axr.set_ylim(0, 1.18)
    axr.set_ylabel('single-gene power', fontsize=7.0)
    axr.tick_params(axis='y', labelsize=6.4)
    axr.spines['top'].set_visible(False)

    # 图例拆两个（2026-09-21）：数据线图例沿用 09-20 版的 lower right 原位不变；
    # 两条参考线单独放左上空白区。合并在一个 6 行图例会太高，实测红线会穿过
    # "BH boundary |Z| at q = 0.05" 这行的文字。
    h1, l1 = axa.get_legend_handles_labels()
    h2, l2 = axr.get_legend_handles_labels()
    d_ = dict(list(zip(l1, h1)) + list(zip(l2, h2)))
    data_order = [L_BH] + power_labs
    ref_order = [L_Z3, L_P8]
    missing = [k for k in data_order + ref_order if k not in d_]
    if missing:
        print('!! 图例句柄缺失:', missing)
        return 1
    lg1 = axa.legend([d_[k] for k in data_order], data_order, frameon=False, fontsize=6.2,
                     loc='lower right', handlelength=1.8, borderpad=0.2)
    axa.add_artist(lg1)
    axa.legend([d_[k] for k in ref_order], ref_order, frameon=False, fontsize=6.0,
               loc='upper left', handlelength=1.8, borderpad=0.2, handletextpad=0.6)
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
