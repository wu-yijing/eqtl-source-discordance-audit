#!/usr/bin/env python3
"""
方案2 可视化脚本: Love plot, Density plot, Scatter plot
"""

import os, json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False

BASE = r"E:\workbuddy\car dia 投稿资料"
RE_DESIGN = os.path.join(BASE, "重新设计方案")
OUTPUT = r"E:\workbuddy\2026-06-27-15-13-34"

# ========== 读取数据 ==========
covar = pd.read_csv(os.path.join(BASE, r"2.候选蛋白协变量s-prediXan\Table_S1_Covariate_Matrix_FINAL_v2.csv"))
matches = pd.read_csv(os.path.join(RE_DESIGN, r"3.富集分析基础数据不完整 + 协变量不匹配\mahalanobis_matched_pairs.csv"))
viz = pd.read_csv(os.path.join(OUTPUT, "viz_z_distribution.csv"))
eqtlgen = pd.read_csv(os.path.join(RE_DESIGN, r"3.富集分析基础数据不完整 + 协变量不匹配\eqtlgen_spredixcan_results.csv"))
comp = pd.read_csv(os.path.join(RE_DESIGN, r"3.富集分析基础数据不完整 + 协变量不匹配\eqtlgen_vs_gtex_comparison.csv"))

# ========== FIG 1: Love Plot — 匹配前后协变量平衡性 ==========
print("  Fig 1: Love plot...")

# 从 matches 提取匹配配对: candidate=1 vs control=0
cand_match = matches[matches['treated'] == 1].copy()
ctrl_match = matches[matches['treated'] == 0].copy()

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 匹配前: 从协变量表中提取
covar_cand = covar[covar['Group'].str.contains('Candidate', case=False)]
covar_non = covar[covar['Group'].str.contains('NonCandidate', case=False)]
covar_t2dm = covar[covar['Group'].str.contains('T2DM', case=False)]

# 只取匹配的候选和对应的对照
cand_genes = set(cand_match['Gene'])
ctrl_genes = set(ctrl_match['Gene'])
covar_cand_matched = covar[covar['Gene'].isin(cand_genes)]
covar_ctrl_matched = covar[covar['Gene'].isin(ctrl_genes)]

# 标准化均值差异 (SMD)
def calc_smd(df1, df2, col):
    m1, m2 = df1[col].mean(), df2[col].mean()
    s1, s2 = df1[col].std(), df2[col].std()
    n1, n2 = len(df1), len(df2)
    # pooled sd
    sp = np.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2))
    return (m1 - m2) / sp if sp > 0 else 0

# 匹配前 SMD
cols_log = ['Length_bp', 'GC_pct']
cols_eqtl = ['Max_eQTL_SNPs']

smd_before = {}
smd_after = {}

for col in cols_log:
    # 匹配前: 候选 vs 匹配对照
    smd_before[col] = calc_smd(covar_cand, covar_ctrl_matched, col)
    # 匹配后: 通过匹配对计算
    merged = pd.merge(cand_match, ctrl_match, on='subclass', suffixes=('_cand', '_ctrl'))
    diff = merged[col.replace('_bp','') + '_cand'] - merged[col.replace('_bp','') + '_ctrl'] if col != 'GC_pct' else \
           merged['Length_bp_cand'] - merged['Length_bp_ctrl']
    # 转换为 SMD: 匹配后的差值应接近 0
    mean_diff = merged[f'log10_Length_cand'].mean() - merged[f'log10_Length_ctrl'].mean() if col == 'Length_bp' else \
                merged['Length_bp_cand'].mean() - merged['Length_bp_ctrl'].mean()
    smd_after[col] = mean_diff / covar[col].std() if covar[col].std() > 0 else 0

# 使用 match 中的 log10_Length 计算
merged_len = pd.merge(cand_match, ctrl_match, on='subclass', suffixes=('_cand', '_ctrl'))
# log10 Length
diff_loglen = merged_len['log10_Length_cand'] - merged_len['log10_Length_ctrl']
smd_after['Length_bp'] = diff_loglen.mean() / merged_len['log10_Length_cand'].std() if merged_len['log10_Length_cand'].std() > 0 else 0

# n_eQTL_SNPs
diff_snp = merged_len['n_eQTL_SNPs_cand'] - merged_len['n_eQTL_SNPs_ctrl']
smd_after['n_eQTL_SNPs'] = diff_snp.mean() / merged_len['n_eQTL_SNPs_cand'].std() if merged_len['n_eQTL_SNPs_cand'].std() > 0 else 0

# GC_pct - 用协变量表的
cand_gc = covar_cand.set_index('Gene')['GC_pct']
ctrl_gc = covar_ctrl_matched.set_index('Gene')['GC_pct']
common_gc = cand_gc.index.intersection(ctrl_gc.index)
if len(common_gc) > 0:
    diff_gc = cand_gc.loc[common_gc].values - ctrl_gc.loc[common_gc].values
    smd_after['GC_pct'] = diff_gc.mean() / covar['GC_pct'].std()
else:
    smd_after['GC_pct'] = 0

# Plot
params = ['Length_bp', 'GC_pct', 'n_eQTL_SNPs']
colors = ['#E74C3C', '#3498DB']
for i, (ax, p) in enumerate(zip(axes, params)):
    before = [smd_before.get(p, 0)]
    after = [smd_after.get(p, 0)]
    y = [0]
    ax.scatter(before, y, color=colors[0], s=100, marker='o', label='匹配前', zorder=5)
    ax.scatter(after, y, color=colors[1], s=100, marker='s', label='匹配后', zorder=5)
    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(-0.1, color='gray', linestyle=':', alpha=0.3)
    ax.axvline(0.1, color='gray', linestyle=':', alpha=0.3)
    ax.set_xlabel('标准化均值差 (SMD)')
    ax.set_title(f'{p}')
    ax.set_yticks([])
    # Add text annotation
    ax.text(0.02, 0.9, f'前={before[0]:.3f}\n后={after[0]:.3f}', transform=ax.transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    if i == 0:
        ax.legend(loc='lower right', fontsize=8)

fig.suptitle('Mahalanobis 匹配前后协变量平衡性 (Love Plot)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT, "Fig1_Love_Plot.png"), dpi=200, bbox_inches='tight')
plt.close()
print(f"  → 已保存: Fig1_Love_Plot.png")

# ========== FIG 2: Density plot — 三组 |Z| 分布 ==========
print("  Fig 2: Density plot...")

# 使用 eQTLGen 数据
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
        sub = eqtlgen[(eqtlgen['Trait'] == trait) & (eqtlgen['Group'] == grp)]['Abs_Z'].dropna()
        if len(sub) < 2:
            continue
        # 使用核密度估计
        from scipy.stats import gaussian_kde
        try:
            kde = gaussian_kde(sub.values)
            x = np.linspace(0, max(sub.values) + 1, 100)
            ax.plot(x, kde(x), color=group_colors[grp], label=group_labels[grp], linewidth=2)
            ax.fill_between(x, kde(x), alpha=0.1, color=group_colors[grp])
        except:
            ax.hist(sub.values, bins=10, density=True, alpha=0.5, color=group_colors[grp], label=group_labels[grp])
    
    ax.set_xlabel('|Z| score (eQTLGen)')
    ax.set_ylabel('密度')
    ax.set_title(f'{trait} — |Z| 分布')
    ax.legend(fontsize=8)
    ax.axvline(1.96, color='red', linestyle='--', alpha=0.5, label='|Z|=1.96')

fig.suptitle('三组 eQTLGen TWAS |Z| 密度分布', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT, "Fig2_Density_Plot.png"), dpi=200, bbox_inches='tight')
plt.close()
print(f"  → 已保存: Fig2_Density_Plot.png")

# ========== FIG 3: Scatter plot — GTEx vs eQTLGen Z 值对比 ==========
print("  Fig 3: Scatter plot...")

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

scatter_colors = {
    'Candidate': '#E74C3C',
    'NonCandidate': '#3498DB',
    'T2DM_Control': '#2ECC71'
}

for i, (ax, trait) in enumerate(zip(axes, ['DR', 'DN', 'DPN'])):
    sub = comp[comp['Trait'] == trait].copy()
    if len(sub) == 0:
        ax.text(0.5, 0.5, '无数据', ha='center', va='center', fontsize=14)
        continue
    
    # 按分组着色
    for grp, color in scatter_colors.items():
        grp_sub = sub[sub['Group'] == grp]
        if len(grp_sub) > 0:
            same = grp_sub[grp_sub['Same_Direction'] == True]
            diff = grp_sub[grp_sub['Same_Direction'] == False]
            ax.scatter(same['Z_GTEx'], same['Z_eQTLGen'], color=color, marker='o', alpha=0.6, s=40, label=f'{grp.split("_")[1]} (一致)')
            ax.scatter(diff['Z_GTEx'], diff['Z_eQTLGen'], color=color, marker='x', alpha=0.8, s=40)
    
    # 对角线
    lim = max(sub[['Z_eQTLGen', 'Z_GTEx']].max().max(), abs(sub[['Z_eQTLGen', 'Z_GTEx']].min().min()))
    ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.3, linewidth=1)
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.3)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.3)
    
    ax.set_xlabel('Z_GTEx')
    ax.set_ylabel('Z_eQTLGen')
    ax.set_title(f'{trait} — GTEx vs eQTLGen Z 值')
    ax.legend(fontsize=7, loc='upper left')
    ax.set_xlim(-lim*1.1, lim*1.1)
    ax.set_ylim(-lim*1.1, lim*1.1)

fig.suptitle('GTEx vs eQTLGen TWAS Z 值散点图', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT, "Fig3_Scatter_Plot.png"), dpi=200, bbox_inches='tight')
plt.close()
print(f"  → 已保存: Fig3_Scatter_Plot.png")

# ========== FIG 4: 方向一致性柱状图 ==========
print("  Fig 4: Direction consistency bar chart...")

fig, ax = plt.subplots(figsize=(10, 6))

# 按组计算方向一致性
groups_ordered = ['Candidate', 'NonCandidate', 'T2DM_Control']
bar_data = []
for grp in groups_ordered:
    sub = comp[comp['Group'] == grp]
    if len(sub) == 0:
        continue
    n_same = sub['Same_Direction'].sum()
    n_total = len(sub)
    pct = n_same / n_total * 100
    bar_data.append({'Group': group_labels[grp], '一致': n_same, '不一致': n_total - n_same, '一致率': pct})

df_bar = pd.DataFrame(bar_data)
x = np.arange(len(df_bar))
width = 0.35

bars1 = ax.bar(x - width/2, df_bar['一致'], width, label='方向一致', color='#2ECC71')
bars2 = ax.bar(x + width/2, df_bar['不一致'], width, label='方向不一致', color='#E74C3C')

# 添加百分比标签
for i, (_, r) in enumerate(df_bar.iterrows()):
    ax.text(i, r['一致'] + 0.5, f"{r['一致率']:.1f}%", ha='center', fontsize=10, fontweight='bold')

ax.set_xlabel('分组')
ax.set_ylabel('基因数')
ax.set_title('GTEx vs eQTLGen 方向一致性 (所有表型合并)')
ax.set_xticks(x)
ax.set_xticklabels(df_bar['Group'], fontsize=9)
ax.legend()
ax.axhline(len(df_bar) / 2, color='gray', linestyle='--', alpha=0.3, label='随机期望50%')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT, "Fig4_Direction_Consistency.png"), dpi=200, bbox_inches='tight')
plt.close()
print(f"  → 已保存: Fig4_Direction_Consistency.png")

print("\n所有可视化完成!")
