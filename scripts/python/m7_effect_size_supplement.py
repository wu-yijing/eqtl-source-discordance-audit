#!/usr/bin/env python3
"""
M7 Effect Size Supplement
==========================
Compute:
1. Top-k overlap analysis (GTEx vs eQTLGen, by |Z| ranking)
2. Fisher's z 95% CI for decomposition arm Spearman ρ values
3. Cross-source Z-score concordance by phenotype

Output: JSON for manuscript integration
"""
import sys
sys.path.insert(0, r'C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Lib/site-packages')
import pandas as pd
import numpy as np
from scipy import stats
from math import sqrt, atanh, tanh
import json, os

REPO = r"E:/workbuddy/eqtl-source-discordance-audit"
DATA = os.path.join(REPO, "data", "processed")
# ⚠️ SCZ 快照属 in-house 估计器（n = 2,511），已被官方 MetaXcan v0.8.1 取代，**不随本仓库分发**。
#    如只需重跑主臂（primary arm），本项可缺省；需要 SCZ 部分时用 M7_SCZ_SNAPSHOT 指向你的副本。
SCZ = os.environ.get("M7_SCZ_SNAPSHOT", os.path.join(REPO, "scz_results"))

# ============================================================
# 1. Load primary comparison data
# ============================================================
comp = pd.read_csv(os.path.join(DATA, "eqtlgen_vs_gtex_comparison.csv"))
print(f"Primary comparison: {len(comp)} rows, columns: {list(comp.columns)}")

# ============================================================
# 2. Top-k overlap analysis
# ============================================================
def top_k_overlap(df, z1_col, z2_col, ks=[5, 10, 20, 30]):
    """Compute overlap and Jaccard index for top-k genes by |Z|."""
    results = []
    valid = df[[z1_col, z2_col]].dropna()
    n = len(valid)
    
    rank1 = valid[z1_col].abs().sort_values(ascending=False)
    rank2 = valid[z2_col].abs().sort_values(ascending=False)
    
    for k in ks:
        top1 = set(rank1.head(k).index)
        top2 = set(rank2.head(k).index)
        overlap = len(top1 & top2)
        union = len(top1 | top2)
        jaccard = overlap / union if union > 0 else 0
        results.append({
            'k': k,
            'n_total': n,
            'overlap': overlap,
            'jaccard': round(jaccard, 3),
            'overlap_pct': round(100 * overlap / k, 1),
            'gt_gain': len(top2 - top1),
            'eqtlgen_gain': len(top1 - top2)
        })
    return results

# Overall top-k overlap
overall_topk = top_k_overlap(comp, 'Z_GTEx', 'Z_eQTLGen')
print("\n=== Top-k overlap (overall, n=102) ===")
for r in overall_topk:
    print(f"  k={r['k']:2d}: overlap={r['overlap']}, Jaccard={r['jaccard']}, %={r['overlap_pct']}%")

# Top-k overlap by phenotype
pheno_topk = {}
for pheno in comp['Trait'].unique():
    df_p = comp[comp['Trait'] == pheno]
    pheno_topk[pheno] = top_k_overlap(df_p, 'Z_GTEx', 'Z_eQTLGen')
    print(f"\n  {pheno} (n={len(df_p)}):")
    for r in pheno_topk[pheno]:
        print(f"    k={r['k']:2d}: overlap={r['overlap']}, Jaccard={r['jaccard']}")

# ============================================================
# 3. Fisher's z 95% CI for Spearman ρ
# ============================================================
def fisher_z_ci(rho, n, alpha=0.05):
    """Compute Fisher's z 95% CI for Spearman/Bravais-Pearson correlation."""
    z = atanh(rho)
    se = 1.0 / sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha/2)
    lo = tanh(z - z_crit * se)
    hi = tanh(z + z_crit * se)
    return round(lo, 3), round(hi, 3)

# Load SCZ decomposition data
with open(os.path.join(SCZ, "scz_decomp_limit0.json")) as f:
    scz_decomp = json.load(f)
with open(os.path.join(SCZ, "scz_axis_difference.json")) as f:
    scz_axis = json.load(f)

# Build Fisher's z CI table for all decomposition arms
print("\n=== Fisher's z 95% CI for decomposition arms ===")
arms = [
    ("Primary (GTEx v8 multi vs eQTLGen)", 0.29, 102, "HOTAIR testbed"),
    ("Panel-only (GTEx WB vs eQTLGen WB)", 0.31, 162, "HOTAIR testbed"),
    ("Tissue-only (GTEx WB vs NT)", 0.45, 150, "HOTAIR testbed"),
    ("Dual mismatch (GTEx multi vs eQTLGen)", 0.26, 201, "HOTAIR testbed"),
    ("SCZ Source axis (eQTLGen vs GTEx multi)", 0.522, 2511, "SCZ replication"),
    ("SCZ Tissue axis (GTEx WB vs NT)", 0.509, 2511, "SCZ replication"),
]

fisher_results = []
for label, rho, n, source in arms:
    lo, hi = fisher_z_ci(rho, n)
    fisher_results.append({
        'arm': label,
        'n': n,
        'spearman_rho': rho,
        'fisher_z_95ci_lo': lo,
        'fisher_z_95ci_hi': hi,
        'source': source
    })
    print(f"  {label}: ρ={rho}, n={n}, 95% CI = [{lo}, {hi}]")

# ============================================================
# 4. Direction consistency by Z-score magnitude strata
# ============================================================
print("\n=== Direction consistency by Z-score threshold ===")
def dir_consistency_by_z(df, z1_col, z2_col, thresholds=[1.0, 1.96, 2.58, 3.0]):
    valid = df[[z1_col, z2_col, 'Same_Direction']].dropna()
    results = []
    
    for thresh in thresholds:
        # Both |Z| >= threshold
        both_strong = valid[(valid[z1_col].abs() >= thresh) & (valid[z2_col].abs() >= thresh)]
        both_dir = both_strong['Same_Direction'].sum()
        both_n = len(both_strong)
        
        # Only GTEx |Z| >= threshold
        gtex_strong = valid[valid[z1_col].abs() >= thresh]
        gtex_dir = gtex_strong['Same_Direction'].sum()
        gtex_n = len(gtex_strong)
        
        # Only eQTLGen |Z| >= threshold
        eqtl_strong = valid[valid[z2_col].abs() >= thresh]
        eqtl_dir = eqtl_strong['Same_Direction'].sum()
        eqtl_n = len(eqtl_strong)
        
        results.append({
            'threshold': thresh,
            'both_n': both_n,
            'both_dir_pct': round(100*both_dir/both_n, 1) if both_n > 0 else None,
            'gtex_strong_n': gtex_n,
            'gtex_strong_dir_pct': round(100*gtex_dir/gtex_n, 1) if gtex_n > 0 else None,
            'eqtlgen_strong_n': eqtl_n,
            'eqtlgen_strong_dir_pct': round(100*eqtl_dir/eqtl_n, 1) if eqtl_n > 0 else None
        })
    return results

dir_strata = dir_consistency_by_z(comp, 'Z_GTEx', 'Z_eQTLGen')
for r in dir_strata:
    print(f"  |Z|≥{r['threshold']}: both_n={r['both_n']}, both_dir={r['both_dir_pct']}%; "
          f"GTEx_strong_dir={r['gtex_strong_dir_pct']}%; eQTLGen_strong_dir={r['eqtlgen_strong_dir_pct']}%")

# ============================================================
# 5. Save all results
# ============================================================
output = {
    'top_k_overlap_overall': overall_topk,
    'top_k_overlap_by_phenotype': pheno_topk,
    'fisher_z_ci': fisher_results,
    'direction_consistency_by_z_strata': dir_strata,
    'note': 'Computed 2026-08-09 from data/processed/eqtlgen_vs_gtex_comparison.csv (primary, n=102) and the superseded in-house SCZ snapshot (not distributed with this repository)'
}

outpath = os.path.join(SCZ, "..", "m7_effect_size_supplement.json")
with open(outpath, 'w') as f:
    json.dump(output, f, indent=2, ensure_ascii=False)
print(f"\nResults saved to: {outpath}")
