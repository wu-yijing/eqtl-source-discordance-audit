#!/usr/bin/env python3
"""
Regenerate Figure 2 with legend moved below the panels and above the title.
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'processed')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'figs')
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'font.size': 11,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.unicode_minus': False,
})

C_NONCAND = '#3498DB'
C_T2DM = '#2ECC71'
C_CANDIDATE = '#E74C3C'

def load_comp():
    comp = pd.read_csv(os.path.join(DATA_DIR, 'eqtlgen_vs_gtex_comparison.csv'))
    comp['Group'] = comp['Group'].replace({
        '30_HOTAIR_Candidate': 'Candidate',
        '44_NonCandidate_HOTAIR': 'NonCandidate',
        '30_T2DM_Control': 'T2DM_Control',
        'Non-Candidate': 'NonCandidate',
        'T2DM Control': 'T2DM_Control',
    })
    return comp

def generate_figure2():
    comp = load_comp()
    comp = comp[comp['Z_eQTLGen'].notna() & comp['Z_GTEx'].notna()]

    fig, axes = plt.subplots(1, 3, figsize=(6.69, 2.90))
    fig.subplots_adjust(bottom=0.40, top=0.88, wspace=0.38)

    group_colors = {
        'Candidate': C_CANDIDATE,
        'NonCandidate': C_NONCAND,
        'T2DM_Control': C_T2DM,
    }
    group_labels = {
        'Candidate': 'Candidates',
        'NonCandidate': 'Non-Candidates',
        'T2DM_Control': 'T2DM Controls',
    }

    for i, (ax, trait) in enumerate(zip(axes, ['DR', 'DN', 'DPN'])):
        sub = comp[comp['Trait'] == trait]
        if len(sub) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                    transform=ax.transAxes, fontsize=12, color='gray')
            continue

        rho, pval = spearmanr(sub['Z_GTEx'], sub['Z_eQTLGen'])
        rho_str = f'ρ = {rho:.3f}' + ('*' if pval < 0.001 else '')

        for grp, color in group_colors.items():
            grp_sub = sub[sub['Group'] == grp]
            if len(grp_sub) == 0:
                continue
            same = grp_sub[grp_sub['Same_Direction'] == True]
            diff = grp_sub[grp_sub['Same_Direction'] != True]

            ax.scatter(same['Z_GTEx'], same['Z_eQTLGen'],
                       c=color, marker='o', alpha=0.5, s=28,
                       edgecolors='none')
            ax.scatter(diff['Z_GTEx'], diff['Z_eQTLGen'],
                       c=color, marker='x', alpha=0.7, s=32, linewidths=1.0)

        lim = max(abs(sub[['Z_eQTLGen', 'Z_GTEx']].values.flatten()).max(), 1.5)
        ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.25, lw=0.8)
        ax.axhline(0, color='gray', lw=0.4, alpha=0.3)
        ax.axvline(0, color='gray', lw=0.4, alpha=0.3)

        ax.set_xlabel('')  # shared xlabel below panels
        ax.set_ylabel('Z-score (eQTLGen)', fontsize=10)
        ax.set_title(f'{trait}', fontsize=12, fontweight='bold')
        ax.set_xlim(-lim * 1.1, lim * 1.1)
        ax.set_ylim(-lim * 1.1, lim * 1.1)

        count_text = f'n = {len(sub)}'
        ax.text(0.97, 0.03, f'{count_text}\n{rho_str}',
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=8, fontstyle='italic',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.5))

    # Shared x-axis label below panels
    fig.text(0.5, 0.30, 'Z-score (GTEx v8)', fontsize=10, ha='center', va='center')

    # Shared legend below the x-label, centered
    legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=C_NONCAND,
                   markersize=7, alpha=0.5, label='Non-Candidates (same)',
                   markeredgecolor='none', linestyle='None'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=C_T2DM,
                   markersize=7, alpha=0.5, label='T2DM Controls (same)',
                   markeredgecolor='none', linestyle='None'),
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=2,
               frameon=True, framealpha=0.95, fontsize=8,
               bbox_to_anchor=(0.5, 0.19), columnspacing=1.2, handletextpad=0.3)

    # Title below the legend, smaller font
    fig.text(0.5, 0.08, 'Figure 2. GTEx vs eQTLGen TWAS Z-score Comparison',
             fontsize=11, fontweight='bold', ha='center', va='bottom')

    out_pdf = os.path.join(FIG_DIR, 'Figure2.pdf')
    out_png = os.path.join(FIG_DIR, 'Figure2.png')
    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close()
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")

if __name__ == '__main__':
    generate_figure2()
