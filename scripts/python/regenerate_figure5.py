#!/usr/bin/env python3
"""
Regenerate Figure 5 (RNH1 effect attenuation) with:
- Main title moved below the figure, centered, smaller font
- Subtitles (a) and (b) placed below each panel, centered, smaller font
- Panel (a) legend shrunk and placed in the upper-right corner
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'processed')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'figs')
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'font.size': 8,
    'axes.titlesize': 8,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 6,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.unicode_minus': False,
})

C_CANDIDATE = '#E74C3C'
C_NONCAND = '#3498DB'
C_CANDIDATE_LIGHT = '#FADBD8'
C_NONCAND_LIGHT = '#D6EAF8'


def load_data():
    eqtlgen = pd.read_csv(os.path.join(DATA_DIR, 'eqtlgen_spredixcan_results.csv'))
    eqtlgen['Group'] = eqtlgen['Group'].replace({
        '30_HOTAIR_Candidate': 'Candidate',
        '44_NonCandidate_HOTAIR': 'NonCandidate',
        '30_T2DM_Control': 'T2DM_Control',
    })
    rnh1_gtex = {
        'RNH1': {'DR': 13.82, 'DN': 7.00, 'DPN': 8.57}
    }
    return eqtlgen, rnh1_gtex


def generate_figure5():
    eqtlgen, rnh1_gtex = load_data()

    fig, axes = plt.subplots(1, 2, figsize=(6.69, 2.95))
    fig.subplots_adjust(bottom=0.30, top=0.82, left=0.10, right=0.95, wspace=0.30)

    # --- Panel A: GTEx vs eQTLGen ---
    ax = axes[0]
    traits = ['DR', 'DN', 'DPN']
    x = np.arange(len(traits))
    width = 0.30

    gtex_vals = [rnh1_gtex['RNH1'][t] for t in traits]
    eqtlgen_vals = []
    for t in traits:
        v = eqtlgen[(eqtlgen['Gene'] == 'RNH1') & (eqtlgen['Trait'] == t)]['Z_eQTLGen'].values
        eqtlgen_vals.append(abs(v[0]) if len(v) > 0 else 0)

    bars1 = ax.bar(x - width/2, gtex_vals, width, label='GTEx v8',
                   color=C_CANDIDATE_LIGHT, edgecolor=C_CANDIDATE, linewidth=0.8)
    bars2 = ax.bar(x + width/2, eqtlgen_vals, width, label='eQTLGen',
                   color=C_NONCAND_LIGHT, edgecolor=C_NONCAND, linewidth=0.8)

    for bars, vals in [(bars1, gtex_vals), (bars2, eqtlgen_vals)]:
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                    f'|Z|={val:.1f}', ha='center', fontsize=6, fontweight='bold')

    ax.axhline(1.96, color='red', linestyle='--', alpha=0.4, lw=0.8)
    ax.text(2.4, 2.0, '|Z|=1.96', fontsize=7, color='red')
    ax.set_xticks(x)
    ax.set_xticklabels(traits, fontsize=7)
    ax.set_ylabel('|Z-score|', fontsize=8)
    ax.set_ylim(0, max(max(gtex_vals), max(eqtlgen_vals)) * 1.3)
    ax.legend(fontsize=6, loc='upper right', framealpha=0.9,
              bbox_to_anchor=(0.98, 0.98), borderpad=0.2, labelspacing=0.15,
              handletextpad=0.25, handlelength=1.0)

    # --- Panel B: Cross-population replication ---
    ax = axes[1]
    cohorts = ['FinnGen\nR13 DR', 'UKB\nGCST90043640', 'UKB\nXue 2022',
               'FinnGen\nDN', 'Sakaue\n2021 DN']
    z_vals = [8.50, 0.57, 0.95, 4.54, 3.59]
    ci_low = [8.26, 0.33, 0.68, 4.30, 3.20]
    ci_high = [8.74, 0.81, 1.22, 4.78, 3.98]
    colors = [C_CANDIDATE, C_NONCAND, C_NONCAND, C_CANDIDATE, C_NONCAND]

    y = np.arange(len(cohorts))
    for i, (cohort, z, cl, ch, col) in enumerate(zip(cohorts, z_vals, ci_low, ci_high, colors)):
        ax.errorbar(z, i, xerr=[[z-cl], [ch-z]], fmt='o', color=col,
                    markersize=6, capsize=3, capthick=1.2, linewidth=1.2)
        ax.text(z + 0.3, i, f'Z={z:.2f}', va='center', fontsize=6,
                fontweight='bold' if z > 1.96 else 'normal')

    ax.axvline(0, color='gray', linestyle='-', alpha=0.3, lw=0.5)
    ax.axvline(1.96, color='red', linestyle='--', alpha=0.4, lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(cohorts, fontsize=6)
    ax.set_xlabel('Z-score', fontsize=8)
    ax.set_xlim(-1, 10)
    ax.invert_yaxis()

    # --- Subtitles below each panel, centered ---
    fig.text(0.28, 0.16, '(a) RNH1: GTEx vs eQTLGen',
             fontsize=9, fontweight='bold', ha='center', va='center')
    fig.text(0.78, 0.16, '(b) RNH1 Cross-Population Replication',
             fontsize=9, fontweight='bold', ha='center', va='center')

    # --- Main title below the figure, centered ---
    fig.text(0.5, 0.03,
             'Figure 5. RNH1 Effect Attenuation Across eQTL Sources and Populations',
             fontsize=10, fontweight='bold', ha='center', va='bottom')

    plt.savefig(os.path.join(FIG_DIR, 'Figure5.png'))
    plt.savefig(os.path.join(FIG_DIR, 'Figure5.pdf'))
    plt.close()
    print('Figure 5 saved to', FIG_DIR)


if __name__ == '__main__':
    generate_figure5()
