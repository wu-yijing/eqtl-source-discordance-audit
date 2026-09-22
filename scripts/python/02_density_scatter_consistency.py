#!/usr/bin/env python3
"""
方案2 可视化脚本 (v2 - 简化版)
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False

# Use repository-relative paths for portability
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'processed')
OUTPUT = os.path.join(REPO_ROOT, 'analyses')

# 读取数据 — use repo data/processed/ files where available
covar = pd.read_csv(os.path.join(DATA_DIR, "covariate_matrix.csv"))
matches = pd.read_csv(os.path.join(DATA_DIR, "mahalanobis_matched_pairs.csv"))
comp = pd.read_csv(os.path.join(DATA_DIR, "eqtlgen_vs_gtex_comparison.csv"))
eqtlgen = pd.read_csv(os.path.join(DATA_DIR, "eqtlgen_spredixcan_results.csv"))

# ====== Fig 1: Love Plot (简化版) ======
print("  Fig 1: Love Plot...")
# 匹配前的候选基因
cand_genes = covar[covar['Group'].str.contains('Candidate', case=False)]
# 匹配后的对照基因
ctrl_genes = set(matches[matches['treated'] == 0]['Gene'])
ctrl_covar = covar[covar['Gene'].isin(ctrl_genes)]

fig, axes = plt.subplots(1, 3, figsize=(14, 5))
params = [('Length_bp', '基因长度 (bp, log10)', True),
          ('GC_pct', 'GC含量 (%)', False),
          ('Max_eQTL_SNPs', 'eQTL SNP数', False)]
colors_before = '#E74C3C'
colors_after = '#3498DB'

for i, (col, label, log_scale) in enumerate(params):
    # 匹配前: 候选 vs 匹配对照
    before_cand = cand_genes[col].dropna()
    before_ctrl = ctrl_covar[col].dropna()
    if log_scale:
        before_cand = np.log10(before_cand.replace(0, np.nan)).dropna()
        before_ctrl = np.log10(before_ctrl.replace(0, np.nan)).dropna()
    
    # 匹配后: 从matches取配对
    cand_match = matches[matches['treated'] == 1].set_index('subclass')
    ctrl_match = matches[matches['treated'] == 0].set_index('subclass')
    merged = cand_match.join(ctrl_match, lsuffix='_cand', rsuffix='_ctrl')
    
    if col == 'Length_bp':
        after_cand = merged['log10_Length_cand']
        after_ctrl = merged['log10_Length_ctrl']
    elif col == 'Max_eQTL_SNPs':
        after_cand = merged['n_eQTL_SNPs_cand']
        after_ctrl = merged['n_eQTL_SNPs_ctrl']
    else:
        # GC_pct - 从协变量表取
        gc_dict = dict(zip(covar['Gene'], covar['GC_pct']))
        after_cand = merged['Gene_cand'].map(gc_dict)
        after_ctrl = merged['Gene_ctrl'].map(gc_dict)
    
    # 绘制匹配前分布
    bp1 = ax_box = axes[i]
    bp1.boxplot([before_cand, before_ctrl], positions=[0, 1], widths=0.5,
                patch_artist=True,
                boxprops=dict(facecolor=colors_before, alpha=0.3),
                medianprops=dict(color='red'))
    
    # 绘制匹配后分布
    bp2 = axes[i].boxplot([after_cand.dropna(), after_ctrl.dropna()], positions=[2, 3], widths=0.5,
                          patch_artist=True,
                          boxprops=dict(facecolor=colors_after, alpha=0.3),
                          medianprops=dict(color='red'))
    
    axes[i].set_xticks([0.5, 2.5])
    axes[i].set_xticklabels(['匹配前', '匹配后'], fontsize=9)
    axes[i].set_ylabel(label, fontsize=9)
    axes[i].set_title(f'{col}', fontsize=11)
    
    # 添加样本数
    ylim_bottom, ylim_top = axes[i].get_ylim()
    axes[i].text(0, ylim_bottom, f'n={len(before_cand)}', ha='center', fontsize=8)
    axes[i].text(3, ylim_bottom, f'n={len(after_cand.dropna())}', ha='center', fontsize=8)

fig.suptitle('Mahalanobis 匹配前后协变量平衡性\n(红=候选组, 蓝=匹配对照)', fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT, "Fig1_Love_Plot.png"), dpi=200, bbox_inches='tight')
plt.close()
print("  → 已保存: Fig1_Love_Plot.png")

# ====== Fig 2: Density Plot ======
print("  Fig 2: Density Plot...")
from scipy.stats import gaussian_kde

eqtlgen['Group'] = eqtlgen['Group'].fillna('Unknown')
eqtlgen['Abs_Z'] = eqtlgen['Z_eQTLGen'].abs()

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
group_colors = {
    'Candidate': '#E74C3C',
    'NonCandidate': '#3498DB',
    'T2DM_Control': '#2ECC71'
}
group_labels = {
    'Candidate': 'Candidates',
    'NonCandidate': 'Non-Candidates',
    'T2DM_Control': 'T2DM Controls'
}

for i, (ax, trait) in enumerate(zip(axes, ['DR', 'DN', 'DPN'])):
    for grp in ['Candidate', 'NonCandidate', 'T2DM_Control']:
        sub = eqtlgen[(eqtlgen['Trait'] == trait) & (eqtlgen['Group'] == grp)]['Abs_Z'].dropna().values
        if len(sub) < 2:
            continue
        try:
            kde = gaussian_kde(sub)
            x = np.linspace(max(0, sub.min()-0.5), sub.max() + 1, 100)
            ax.plot(x, kde(x), color=group_colors[grp], label=group_labels[grp], linewidth=2)
            ax.fill_between(x, kde(x), alpha=0.1, color=group_colors[grp])
        except:
            ax.hist(sub, bins=10, density=True, alpha=0.5, color=group_colors[grp], label=group_labels[grp])
    
    ax.set_xlabel('|Z| (eQTLGen)', fontsize=10)
    ax.set_ylabel('密度', fontsize=10)
    ax.set_title(f'{trait}', fontsize=12)
    ax.legend(fontsize=8)
    ax.axvline(1.96, color='red', linestyle='--', alpha=0.5)
    ax.text(2.0, ax.get_ylim()[1]*0.95, '|Z|=1.96', fontsize=8, color='red')

fig.suptitle('三组 eQTLGen TWAS |Z| 密度分布 (eQTLGen Whole Blood)', fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT, "Fig2_Density_Plot.png"), dpi=200, bbox_inches='tight')
plt.close()
print("  → 已保存: Fig2_Density_Plot.png")

# ====== Fig 3: Scatter Plot ======
print("  Fig 3: Scatter Plot...")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

scatter_colors = {
    'Candidate': '#E74C3C',
    'NonCandidate': '#3498DB',
    'T2DM_Control': '#2ECC71'
}
marker_map = {True: 'o', False: 'x'}
alpha_map = {True: 0.6, False: 0.9}

for i, (ax, trait) in enumerate(zip(axes, ['DR', 'DN', 'DPN'])):
    sub = comp[comp['Trait'] == trait].copy()
    if len(sub) == 0:
        ax.text(0.5, 0.5, 'No Data', ha='center', va='center', fontsize=14, transform=ax.transAxes)
        continue
    
    for grp, color in scatter_colors.items():
        grp_sub = sub[sub['Group'] == grp]
        if len(grp_sub) == 0:
            continue
        same = grp_sub[grp_sub['Same_Direction'] == True]
        diff = grp_sub[grp_sub['Same_Direction'] == False]
        ax.scatter(same['Z_GTEx'], same['Z_eQTLGen'], c=color, marker='o', alpha=0.6, s=40, 
                   label=f'{grp.split("_")[1]} (一致)' if trait == 'DR' else '')
        ax.scatter(diff['Z_GTEx'], diff['Z_eQTLGen'], c=color, marker='x', alpha=0.8, s=40)
    
    lim = max(abs(sub[['Z_eQTLGen', 'Z_GTEx']].values.flatten()).max(), 1)
    ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.3, lw=1)
    ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
    ax.axvline(0, color='gray', lw=0.5, alpha=0.3)
    ax.set_xlabel('Z_GTEx', fontsize=10)
    ax.set_ylabel('Z_eQTLGen', fontsize=10)
    ax.set_title(f'{trait}', fontsize=12)
    ax.set_xlim(-lim*1.1, lim*1.1)
    ax.set_ylim(-lim*1.1, lim*1.1)
    if i == 0:
        ax.legend(fontsize=7, loc='upper left')

fig.suptitle('GTEx vs eQTLGen TWAS Z 值对比\n(o=方向一致, x=方向反转)', fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT, "Fig3_Scatter_Plot.png"), dpi=200, bbox_inches='tight')
plt.close()
print("  → 已保存: Fig3_Scatter_Plot.png")

# ====== Fig 4: 方向一致性 ======
print("  Fig 4: Direction Consistency...")
fig, ax = plt.subplots(figsize=(10, 6))

group_labels_short = {'Candidate': 'Candidates', 'NonCandidate': 'Non-Candidates', 'T2DM_Control': 'T2DM Controls'}
bar_data = []
for grp in ['Candidate', 'NonCandidate', 'T2DM_Control']:
    sub = comp[comp['Group'] == grp]
    if len(sub) == 0:
        continue
    n_same = sub['Same_Direction'].sum()
    n_total = len(sub)
    bar_data.append({'Group': group_labels_short[grp], 'Same': n_same, 'Diff': n_total - n_same, 'Pct': n_same/n_total*100})

df_bar = pd.DataFrame(bar_data)
x = np.arange(len(df_bar))
w = 0.35
ax.bar(x - w/2, df_bar['Same'], w, label='方向一致', color='#2ECC71')
ax.bar(x + w/2, df_bar['Diff'], w, label='方向不一致', color='#E74C3C')
for _, r in df_bar.iterrows():
    idx = df_bar.index.get_loc(_)
    ax.text(idx, r['Same'] + 0.5, f"{r['Pct']:.1f}%", ha='center', fontsize=11, fontweight='bold')
ax.set_xlabel('分组', fontsize=11)
ax.set_ylabel('基因数', fontsize=11)
ax.set_title('GTEx vs eQTLGen 方向一致性 (三表型合并)', fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels(df_bar['Group'], fontsize=10)
ax.legend(fontsize=10)
ax.axhline(np.mean(df_bar['Same'] + df_bar['Diff']) * 0.5, color='gray', ls='--', alpha=0.4, label='随机期望50%')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT, "Fig4_Direction_Consistency.png"), dpi=200, bbox_inches='tight')
plt.close()
print("  → 已保存: Fig4_Direction_Consistency.png")

print("\n✅ 所有可视化完成! 4张图已生成。")
