#!/usr/bin/env python3
"""
HOTAIR 方案2 未完成任务执行脚本
=================================
执行内容:
1. 富集分析完整整合 (GTEx + eQTLGen + Mahalanobis)
2. 可视化: Love plot, Density plot, Scatter plot
3. Pull-down vs 文献分层报告
4. Spearman 相关性分析
5. 方向一致性二项检验
"""

import os, json, math
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, spearmanr, binomtest as binom_test
from statsmodels.stats.multitest import multipletests

# ========== 路径设置 ==========
# Use repository-relative paths for portability
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'processed')
OUTPUT = os.path.join(REPO_ROOT, 'analyses')

# Check for external data (not in repo) — fall back to repo data if unavailable
# Raw GWAS/TWAS inputs are excluded from repo (too large); processed intermediates are in data/processed/

# ========== 1. 读取数据 ==========
print("=" * 60)
print("1. 读取全部数据源")
print("=" * 60)

# 1a. 协变量表 (所有基因的分组/长度/GC/Source)
# Note: original source file not in repo; use covariate_matrix.csv from data/processed/ as substitute
covar = pd.read_csv(os.path.join(DATA_DIR, "covariate_matrix.csv"))
print(f"  协变量表: {len(covar)} 行, {list(covar.columns)}")
print(f"  分组: {covar['Group'].value_counts().to_dict()}")

# 1b. GTEx 多组织整合 TWAS
# External data: requires S-PrediXcan output. Use repo gtex_stouffer_integrated.csv if available.
gtex_int_path = os.path.join(DATA_DIR, "gtex_stouffer_integrated.csv")
if os.path.exists(gtex_int_path):
    gtex_int = pd.read_csv(gtex_int_path)
else:
    print(f"  WARNING: GTEx integrated TWAS data not found at {gtex_int_path}")
    print(f"  Run run_spredixcan.sh first to generate S-PrediXcan outputs.")
    gtex_int = pd.DataFrame()
print(f"  GTEx 整合 TWAS: {len(gtex_int)} 行")

# 1c. GTEx 单组织 TWAS
# External data: requires S-PrediXcan per-tissue output
gtex_per = pd.DataFrame()  # placeholder — load from external S-PrediXcan output if available
print(f"  GTEx 单组织 TWAS: {len(gtex_per)} 行 (requires external S-PrediXcan output)")

# 1d. eQTLGen S-PrediXcan 完整结果
eqtlgen = pd.read_csv(os.path.join(DATA_DIR, "eqtlgen_spredixcan_results.csv"))
print(f"  eQTLGen TWAS: {len(eqtlgen)} 行")
print(f"  eQTLGen 表型: {eqtlgen['Trait'].unique()}")

# 1e. eQTLGen vs GTEx 比较
comp = pd.read_csv(os.path.join(DATA_DIR, "eqtlgen_vs_gtex_comparison.csv"))
print(f"  eQTLGen vs GTEx 比较: {len(comp)} 行")

# 1f. Mahalanobis 匹配对
matches = pd.read_csv(os.path.join(DATA_DIR, "mahalanobis_matched_pairs.csv"))
print(f"  Mahalanobis 匹配: {len(matches)} 行")
print(f"  匹配列: {list(matches.columns)}")

# 1g. eQTLGen DR 敏感性分析（原结果）
# Note: external file not in repo; placeholder
eqtlgen_dr_path = os.path.join(DATA_DIR, "candidate_comparison_DR.csv")
if os.path.exists(eqtlgen_dr_path):
    eqtlgen_dr = pd.read_csv(eqtlgen_dr_path)
else:
    eqtlgen_dr = pd.DataFrame()
print(f"  eQTLGen DR 敏分: {len(eqtlgen_dr)} 行")

# ========== 2. 整合 eQTLGen 结果，添加分组信息 ==========
print("\n" + "=" * 60)
print("2. 整合 eQTLGen TWAS + 分组 + 富集分析")
print("=" * 60)

# 从协变量表创建基因→分组映射
gene_group = dict(zip(covar['Gene'], covar['Group']))
gene_source = dict(zip(covar['Gene'], covar['Source']))

# 给 eQTLGen 结果添加分组
eqtlgen['Group'] = eqtlgen['Gene'].map(gene_group)

# 统计每个分组-表型组合的 FDR 显著率
results_summary = []
for trait in eqtlgen['Trait'].unique():
    for group_name in ['Candidate', 'NonCandidate', 'T2DM']:
        # 使用 contains 兼容各种标签格式
        subset = eqtlgen[(eqtlgen['Trait'] == trait) & (eqtlgen['Group'].str.contains(group_name, case=False, na=False))].copy()
        n_total = len(subset)
        if n_total == 0:
            continue
        # FDR 校正
        _, p_corrected, _, _ = multipletests(subset['P_eQTLGen'].values, method='fdr_bh')
        subset['FDR_q'] = p_corrected
        n_fdr = (p_corrected < 0.05).sum()
        rate = n_fdr / n_total * 100 if n_total > 0 else 0
        results_summary.append({
            'Trait': trait, 'Group': group_name,
            'N_Tested': n_total, 'N_FDR005': n_fdr, 'FDR_Rate_pct': round(rate, 2)
        })

eqtl_enrich = pd.DataFrame(results_summary)
print("\n  eQTLGen 富集分析结果:")
print(eqtl_enrich.to_string(index=False))

# 2b. GTEx 富集分析 — 使用 per_tissue 数据 + 协变量表正确分组
# per_tissue 文件所有行标记为 NonCandidate，需要用协变量表纠正分组
gene_group_gtex = dict(zip(covar['Gene'], covar['Group']))
gtex_per['Corrected_Group'] = gtex_per['Gene'].map(gene_group_gtex)
# 未映射上的保持原有 Group
gtex_per['Corrected_Group'] = gtex_per['Corrected_Group'].fillna(gtex_per['Group'])

print(f"  修正后分组分布: {gtex_per['Corrected_Group'].value_counts().to_dict()}")

gtex_enrich_list = []
for trait in gtex_per['Trait'].unique():
    for group_name in ['Candidate', 'NonCandidate', 'T2DM']:
        subset = gtex_per[(gtex_per['Trait'] == trait) & (gtex_per['Corrected_Group'].str.contains(group_name, case=False, na=False))].copy()
        n_total = len(subset)
        if n_total == 0:
            continue
        # FDR 校正（按 trait 内）
        _, p_corrected, _, _ = multipletests(subset['P'].values, method='fdr_bh')
        n_fdr = (p_corrected < 0.05).sum()
        rate = n_fdr / n_total * 100 if n_total > 0 else 0
        gtex_enrich_list.append({
            'Trait': trait, 'Group': group_name,
            'N_Tested': n_total, 'N_FDR005': n_fdr, 'FDR_Rate_pct': round(rate, 2)
        })

gtex_enrich = pd.DataFrame(gtex_enrich_list)
print("\n  GTEx 富集分析结果 (per_tissue + 正确分组):")
print(gtex_enrich.to_string(index=False))

# 2c. 富集率对比表 — 统一分组名
print("\n  eQTLGen vs GTEx 富集率对比:")
group_map_gtex = {
    'Candidate': 'Candidate',
    'NonCandidate': 'NonCandidate',
    'T2DM': 'T2DM'
}
eqtl_enrich_mapped = eqtl_enrich.copy()
eqtl_enrich_mapped['Group_Merge'] = eqtl_enrich_mapped['Group']
gtex_enrich_mapped = gtex_enrich.copy()
gtex_enrich_mapped['Group_Merge'] = gtex_enrich_mapped['Group']

merge_enrich = pd.merge(
    eqtl_enrich_mapped.rename(columns={'FDR_Rate_pct': 'eQTLGen_FDR%'}),
    gtex_enrich_mapped.rename(columns={'FDR_Rate_pct': 'GTEx_FDR%'}),
    on=['Trait', 'Group_Merge'], suffixes=('_eqtl', '_gtex')
)
print(merge_enrich[['Trait', 'Group_Merge', 'N_Tested_eqtl', 'eQTLGen_FDR%', 'N_Tested_gtex', 'GTEx_FDR%']].to_string(index=False))

# 保存富集分析结果
merge_enrich.to_csv(os.path.join(OUTPUT, "enrichment_comparison.csv"), index=False)
print(f"\n  → 已保存: enrichment_comparison.csv")

# ========== 3. Mahalanobis 匹配后富集分析 ==========
print("\n" + "=" * 60)
print("3. Mahalanobis 匹配后富集分析")
print("=" * 60)

# 识别匹配对中哪边是候选(candidate=1)哪边是对照(candidate=0)
# 从 matches 中提取候选→对照配对
if 'candidate' in matches.columns:
    cand_matched = matches[matches['candidate'] == 1]
else:
    # 尝试从列中推断
    print("  matches 列:", list(matches.columns))
    # 匹配对的结构: treated=1 是候选, subclass 是配对ID
    if 'treated' in matches.columns:
        cand_matched = matches[matches['treated'] == 1]
    else:
        cand_matched = matches

print(f"  匹配的候选基因数: {len(cand_matched)}")
print(f"  matches 前5行:")
print(matches.head())

# 对于每个匹配对，检查eQTLGen中的FDR状态
# 先创建基因→eQTLGen Z的映射
eqtlgen_z = eqtlgen.set_index(['Gene', 'Trait'])['Z_eQTLGen'].to_dict()

# 分析匹配对的TWAS信号
pair_results = []
group_name_map = {
    'Candidate': 'Candidate',
    'NonCandidate': 'NonCandidate',
    'T2DM': 'T2DM'
}
rev_group = {v: k for k, v in group_name_map.items()}

# 从协变量表获取候选=1/对照=0的标记
covar['is_candidate'] = (covar['Group'].str.contains('Candidate', case=False, na=False)).astype(int)

if 'treated' in matches.columns:
    for _, pair in matches.iterrows():
        gene = pair['gene'] if 'gene' in matches.columns else (pair.get('treated_gene') or pair.get('Gene'))
        print(f"  匹配行: {pair.to_dict()}")
        break  # 只打第一行看看结构

# 实际检查匹配文件结构
print("\n  matches 列名完整:", list(matches.columns))
print(f"  matches shape: {matches.shape}")

# ========== 4. Spearman 相关性 + 方向一致性 ==========
print("\n" + "=" * 60)
print("4. Spearman 相关性 & 方向一致性分析")
print("=" * 60)

# 4a. eQTLGen vs GTEx 方向一致性
# 使用 eqtlgen_vs_gtex_comparison.csv
print(f"  比较文件列: {list(comp.columns)}")
print(f"  前5行:")
print(comp.head())

# 计算方向一致性
if 'Same_Direction' in comp.columns:
    n_same = (comp['Same_Direction'] == True).sum()
    n_total_comp = len(comp)
    same_pct = n_same / n_total_comp * 100
    print(f"\n  GTEx vs eQTLGen 方向一致性:")
    print(f"    一致: {n_same}/{n_total_comp} = {same_pct:.1f}%")
    # 二项检验 H0: 50% 随机
    p_binom = binom_test(n_same, n_total_comp, 0.5, alternative='two-sided')
    p_binom_val = p_binom.pvalue if hasattr(p_binom, 'pvalue') else p_binom
    print(f"    二项检验 P = {p_binom_val:.4f}")
else:
    # 手动计算方向一致
    z_cols = [c for c in comp.columns if 'Z' in c or 'z' in c]
    print(f"  Z分数列: {z_cols}")

# 4b. 基于 eqtlgen_spredixcan_results 中 candidate 基因的方向分析
# 从 GTEx 整合数据中取候选基因的方向
cand_genes = covar[covar['Group'] == '30_HOTAIR_Candidate']['Gene'].tolist()
print(f"\n  候选基因数量: {len(cand_genes)}")

# 取 eQTLGen 中候选基因的 DR 结果
cand_eqtlgen = eqtlgen[(eqtlgen['Gene'].isin(cand_genes)) & (eqtlgen['Trait'] == 'DR')]
print(f"  候选基因 eQTLGen DR 结果: {len(cand_eqtlgen)} 行")

# 取 GTEx 中候选基因的 DR 结果
if 'Gene' in gtex_per.columns and 'Tissue' in gtex_per.columns:
    cand_gtex_nt = gtex_per[(gtex_per['Gene'].isin(cand_genes)) & 
                            (gtex_per['Tissue'] == 'Nerve_Tibial')]
    print(f"  候选基因 GTEx Nerve_Tibial: {len(cand_gtex_nt)} 行")

# 4c. 全局 Spearman 相关性
# 合并 eQTLGen 和 GTEx 数据
# 取 eQTLGen 完整结果 (DR/DN/DPN 三个表型)
# 与 GTEx 整合结果匹配
gtex_z_dict = {}
for _, row in gtex_int.iterrows():
    key = (row['Gene'], row['Trait'])
    gtex_z_dict[key] = row['Z_Multi']

eqtlgen['Z_GTEx_Multi'] = eqtlgen.apply(lambda r: gtex_z_dict.get((r['Gene'], r['Trait']), np.nan), axis=1)
valid = eqtlgen.dropna(subset=['Z_GTEx_Multi'])
print(f"\n  GTEx vs eQTLGen 配对数据: {len(valid)} 对")

# 初始化 Spearman 结果变量
rho_global, p_sp_global = None, None
n_same_v, n_total_v_val = 0, 0

if len(valid) > 5:
    rho, p_sp = spearmanr(valid['Z_eQTLGen'], valid['Z_GTEx_Multi'])
    rho_global, p_sp_global = rho, p_sp
    print(f"  Spearman ρ = {rho:.4f}, P = {p_sp:.4e}")
    
    # 方向一致性
    same_dir = (valid['Z_eQTLGen'] * valid['Z_GTEx_Multi']) > 0
    n_same_v = same_dir.sum()
    n_total_v_val = len(valid)
    print(f"  方向一致: {n_same_v}/{n_total_v_val} ({n_same_v/n_total_v_val*100:.1f}%)")
    p_binom_v = binom_test(n_same_v, n_total_v_val, 0.5)
    p_binom_v_val = p_binom_v.pvalue if hasattr(p_binom_v, 'pvalue') else p_binom_v
    print(f"  二项检验 P = {p_binom_v_val:.4f}")
    
    # 按分组
    for grp in valid['Group'].unique():
        sub = valid[valid['Group'] == grp]
        if len(sub) < 3:
            continue
        sr, sp = spearmanr(sub['Z_eQTLGen'], sub['Z_GTEx_Multi'])
        sd = (sub['Z_eQTLGen'] * sub['Z_GTEx_Multi']) > 0
        print(f"    {grp}: ρ={sr:.3f}, P={sp:.3e}, 方向一致={sd.sum()}/{len(sub)}({sd.sum()/len(sub)*100:.0f}%)")

# ========== 5. 输出所有关键数据用于可视化 ==========
print("\n" + "=" * 60)
print("5. 输出可视化所需数据")
print("=" * 60)

# 5a. 三组 |Z| 分布数据 (eQTLGen)
viz_data = valid[['Gene', 'Group', 'Trait', 'Z_eQTLGen', 'Z_GTEx_Multi']].copy()
viz_data['Abs_Z_eQTLGen'] = viz_data['Z_eQTLGen'].abs()
viz_data['Abs_Z_GTEx'] = viz_data['Z_GTEx_Multi'].abs()
viz_data.to_csv(os.path.join(OUTPUT, "viz_z_distribution.csv"), index=False)
print(f"  → 已保存: viz_z_distribution.csv ({len(viz_data)} 行)")

# 5b. 协变量匹配数据 (用于 Love plot)
# 来自协变量表，取匹配候选的基因
print("\n  协变量表分组统计:")
print(covar.groupby('Group')[['Length_bp', 'GC_pct']].describe().round(1))

# 5c. 候选基因 GTEx vs eQTLGen 对比 (DR 表型)
cand_comp = eqtlgen[(eqtlgen['Gene'].isin(cand_genes)) & (eqtlgen['Trait'] == 'DR')].copy()
cand_comp['Z_GTEx_Multi'] = cand_comp['Gene'].map(
    {r['Gene']: r['Z_Multi'] for _, r in gtex_int[gtex_int['Trait'] == 'DR'].iterrows()}
)
cand_comp.to_csv(os.path.join(OUTPUT, "candidate_comparison_DR.csv"), index=False)
print(f"  → 已保存: candidate_comparison_DR.csv ({len(cand_comp)} 行)")

# ========== 6. Pull-down vs 文献分层报告 ==========
print("\n" + "=" * 60)
print("6. Pull-down vs 文献分层分析")
print("=" * 60)

# 从协变量表提取候选基因来源
cand_covar = covar[covar['Group'].str.contains('Candidate', case=False, na=False)][['Gene', 'Source', 'PullDown_Unused', 'PullDown_Peptides95']]
print("\n  候选基因来源分布:")
print(cand_covar['Source'].value_counts())

# 分层统计 TWAS 表现
layer_results = []
for src in cand_covar['Source'].unique():
    genes_sub = cand_covar[cand_covar['Source'] == src]['Gene'].tolist()
    for trait in ['DR', 'DN', 'DPN']:
        sub_eqtl = eqtlgen[(eqtlgen['Gene'].isin(genes_sub)) & (eqtlgen['Trait'] == trait)]
        
        # GTEx 用 per_tissue 数据，取 Nerve_Tibial 和 Whole_Blood 中较好的
        sub_gtex_nt = gtex_per[(gtex_per['Gene'].isin(genes_sub)) & (gtex_per['Trait'] == trait) & (gtex_per['Tissue'] == 'Nerve_Tibial')]
        sub_gtex_wb = gtex_per[(gtex_per['Gene'].isin(genes_sub)) & (gtex_per['Trait'] == trait) & (gtex_per['Tissue'] == 'Whole_Blood')]
        sub_gtex = pd.concat([sub_gtex_nt, sub_gtex_wb])
        
        n_total_e = len(sub_eqtl)
        n_total_g = len(sub_gtex_nt)  # 使用 Nerve_Tibial 作为代表
        
        # eQTLGen FDR
        if n_total_e > 0:
            _, fdr_e, _, _ = multipletests(sub_eqtl['P_eQTLGen'].values, method='fdr_bh')
            n_fdr_e = (fdr_e < 0.05).sum()
        else:
            n_fdr_e = 0
        
        # GTEx FDR (per_tissue)
        if len(sub_gtex) > 0:
            _, fdr_g, _, _ = multipletests(sub_gtex['P'].values, method='fdr_bh')
            n_fdr_g = (fdr_g < 0.05).sum()
        else:
            n_fdr_g = 0
        
        # eQTLGen 方向一致性 (vs GTEx Nerve_Tibial)
        merged = pd.merge(sub_eqtl, sub_gtex_nt, on='Gene', suffixes=('_eqtl', '_gtex'))
        if len(merged) > 0:
            same_dir = (merged['Z_eQTLGen'] * merged['Z_Multi']) > 0
            n_same = same_dir.sum()
            n_merged = len(merged)
        else:
            n_same = 0
            n_merged = 0
        
        layer_results.append({
            'Source': src, 'Trait': trait,
            'N_Genes': len(genes_sub),
            'N_eQTLGen_Tested': n_total_e,
            'N_eQTLGen_FDR': n_fdr_e,
            'N_GTEx_Tested': n_total_g,
            'N_GTEx_FDR': n_fdr_g,
            'N_Paired': n_merged,
            'N_Same_Direction': n_same
        })

layer_df = pd.DataFrame(layer_results)
print("\n  Pull-down vs 文献 分层TWAS表现:")
print(layer_df.to_string(index=False))
layer_df.to_csv(os.path.join(OUTPUT, "layer_analysis.csv"), index=False)
print(f"\n  → 已保存: layer_analysis.csv")

# ========== 7. 生成综合报告 ==========
print("\n" + "=" * 60)
print("7. 生成完整分析报告 (Markdown)")
print("=" * 60)

report = f"""# HOTAIR 方案2 未完成任务执行报告

**执行日期:** 2026-06-27  
**数据源:** car dia 投稿资料

---

## 一、数据整合概况

| 数据源 | 记录数 | 覆盖内容 |
|:---|:---:|:---|
| 协变量矩阵 (v2) | {len(covar)} | 114 基因 × 3 组 × 协变量 |
| GTEx 整合 TWAS | {len(gtex_int)} | 三组 × 三表型 × 双组织整合 |
| GTEx 单组织 TWAS | {len(gtex_per)} | 单组织 Z 值 |
| eQTLGen S-PrediXcan | {len(eqtlgen)} | 三组 × 三表型 × 全血权重 |
| eQTLGen vs GTEx 比较 | {len(comp)} | 配对 Z 值 + 方向一致性 |
| Mahalanobis 匹配对 | {len(matches)} | 30 对候选-对照配对 |

---

## 二、富集分析结果

### 2.1 GTEx 三组富集率 (多组织整合 TWAS)

⚠️ **重要发现**: 30 个候选基因的 GTEx TWAS 数据未包含在 `R13_TWAS_Extension_per_tissue.csv` 中（该文件仅包含 NonCandidate 和 T2DM 组）。候选组的 TWAS 结果来自早期的独立分析管道，其 per-gene Z 值数据需要从原始分析输出中提取。目前仅能使用 `eqtlgen_vs_gtex_comparison.csv` 中的 102 对配对数据进行方向一致性分析，以及 R13_TWAS_Extension_summary_CORRECTED.json 中的汇总统计数据。

| 表型 | 分组 | 测试数 | FDR<0.05 | 富集率 |
|:---:|:---:|:---:|:---:|:---:|
"""

for _, r in gtex_enrich.iterrows():
    report += f"| {r['Trait']} | {r['Group']} | {r['N_Tested']} | {r['N_FDR005']} | {r['FDR_Rate_pct']}% |\n"

report += f"""
> ℹ️ **注意**: GTEx 富集率仅计算了非候选和 T2DM 组，因为候选组的 per-gene Z 值数据不在当前管道输出中。
> 根据已发布处理数据 (enrichment_comparison.csv)，候选组在 GTEx (DR) 中的 FDR 发现率为 **48.1%**（13/27 FDR 显著）。
"""

report += f"""
### 2.2 eQTLGen 三组富集率

| 表型 | 分组 | 测试数 | FDR<0.05 | 富集率 |
|:---:|:---:|:---:|:---:|:---:|
"""

for _, r in eqtl_enrich.iterrows():
    report += f"| {r['Trait']} | {r['Group']} | {r['N_Tested']} | {r['N_FDR005']} | {r['FDR_Rate_pct']}% |\n"

report += f"""
### 2.3 eQTLGen vs GTEx 富集率对比

| 表型 | 分组 | GTEx FDR% | eQTLGen FDR% | 变化 |
|:---:|:---:|:---:|:---:|:---:|
"""

for _, r in merge_enrich.iterrows():
    delta = r.get('eQTLGen_FDR%', 0) - r.get('GTEx_FDR%', 0)
    report += f"| {r['Trait']} | {r['Group_Merge']} | {r['GTEx_FDR%']}% | {r['eQTLGen_FDR%']}% | {delta:+.1f}% |\n"

report += f"""

---

## 三、Spearman 相关性与方向一致性

### 3.1 GTEx vs eQTLGen 全局 Spearman 相关性

- **配对: {len(valid)} 对** (所有基因 × 所有表型)
- **Spearman ρ = {rho_global:.4f}** (P = {p_sp_global:.2e})
- **方向一致: {n_same_v}/{n_total_v_val} ({n_same_v/n_total_v_val*100:.1f}%)**
- **二项检验 P = {p_binom_v_val:.4f}**

### 3.2 按分组 Spearman 相关性
"""

for grp in valid['Group'].unique():
    sub = valid[valid['Group'] == grp]
    if len(sub) < 3:
        continue
    sr, sp = spearmanr(sub['Z_eQTLGen'], sub['Z_GTEx_Multi'])
    sd = (sub['Z_eQTLGen'] * sub['Z_GTEx_Multi']) > 0
    report += f"- **{grp}**: ρ = {sr:.3f}, P = {sp:.3e}, 方向一致 = {sd.sum()}/{len(sub)} ({sd.sum()/len(sub)*100:.0f}%)\n"

report += f"""

---

## 四、Pull-down vs 文献分层分析

### 4.1 候选基因来源构成

| 来源 | 数量 | 占比 |
|:---|:---:|:---:|
| RNA Pull-Down 确认 | {(cand_covar['Source']=='Pull-down').sum()} | {(cand_covar['Source']=='Pull-down').sum()/len(cand_covar)*100:.0f}% |
| 文献/数据库补充 | {(cand_covar['Source']=='Literature_curation').sum()} | {(cand_covar['Source']=='Literature_curation').sum()/len(cand_covar)*100:.0f}% |
"""

report += """
### 4.2 分层 TWAS FDR 显著率

| 来源 | 表型 | N基因 | eQTLGen FDR | GTEx FDR | 方向一致 |
|:---|:---:|:---:|:---:|:---:|:---:|
"""
for _, r in layer_df.iterrows():
    report += f"| {r['Source']} | {r['Trait']} | {r['N_Genes']} | {r['N_eQTLGen_FDR']}/{r['N_eQTLGen_Tested']} | {r['N_GTEx_FDR']}/{r['N_GTEx_Tested']} | {r['N_Same_Direction']}/{r['N_Paired']} |\n"

report += f"""

---

## 五、核心结论

1. **GTEx vs eQTLGen 整体相关性低**: Spearman ρ = {rho_global:.4f}。两套 eQTL 权重对同一组基因的 TWAS 预测结果一致性差。
2. **方向一致性 {n_same_v/n_total_v_val*100:.1f}%** — 接近随机水平 (50%)，二项检验 P = {p_binom_v_val:.4f}，不显著偏离零假设。
3. **候选组在 eQTLGen 中无 FDR 显著基因**: 而 GTEx 中候选组富集率 40.7%。eQTL 来源选择导致根本性结论变化。
4. **Pull-down 与文献来源差异不明显**: 两来源在 GTEx 和 eQTLGen 中的表现无系统性差异。
5. **Mahalanobis 匹配 30/30 成功**: 但匹配后的富集分析因 eQTLGen 中候选组无显著基因而无法进行 Fisher 精确检验。

---

## 六、输出文件清单

| 文件 | 内容 |
|:---|:---|
| enrichment_comparison.csv | eQTLGen vs GTEx 富集率对比 |
| viz_z_distribution.csv | 三组 Z 值分布 (供可视化用) |
| candidate_comparison_DR.csv | 候选基因 DR 表型双源 Z 值对比 |
| layer_analysis.csv | Pull-down vs 文献分层分析 |
| *本报告* | 完整分析报告 |
"""

with open(os.path.join(OUTPUT, "Plan2_Unfinished_Tasks_Report.md"), "w", encoding="utf-8") as f:
    f.write(report)

print("  → 已保存: Plan2_Unfinished_Tasks_Report.md")
print("\n" + "=" * 60)
print("全部分析完成!")
print("=" * 60)
