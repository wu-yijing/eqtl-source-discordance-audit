"""
Generate Supplementary Tables S1 and S2 from processed data.
"""

import pandas as pd
import os

DATA = r"E:/workbuddy/eqtl-source-discordance-audit/data/processed"
OUT = 'data/processed'

# ─────────────────────────────────────────────
# Table S1: Gene List with Group Assignment
# ─────────────────────────────────────────────
print("Generating Table S1: Gene List...")

covar = pd.read_csv(os.path.join(DATA, 'covariate_matrix.csv'))

# Clean column names
covar = covar.rename(columns={
    'Gene': 'Gene',
    'Group': 'Group',
    'Source': 'Source',
    'PullDown_Unused': 'PD_Unused_Score',
    'PullDown_Cov_pct': 'PD_Sequence_Coverage_Pct',
    'PullDown_Peptides95': 'PD_Peptides_95pct',
    'Length_bp': 'Gene_Length_bp',
    'GC_pct': 'GC_Content_Pct',
    'eQTL_SNPs_Mean_Num': 'Mean_eQTL_SNP_Count',
})

# Check which candidates were analyzed (had eQTLGen weights)
eqtlgen = pd.read_csv(os.path.join(DATA, 'eqtlgen_spredixcan_results.csv'))
analyzed_genes = set(eqtlgen['Gene'].unique())

# Add status column
covar['Analyzed'] = covar['Gene'].apply(lambda g: 'Yes' if g in analyzed_genes else 'No')

# For T2DM_Control group, relabel
covar['Source'] = covar['Source'].replace({
    'GWAS_literature': 'GWAS Literature',
    'Literature_curation': 'Literature Curation',
    'Pull-down': 'RNA Pull-down'
})

# Select and order columns
t1_cols = ['Gene', 'Group', 'Source', 'Analyzed', 
           'PD_Unused_Score', 'PD_Sequence_Coverage_Pct', 'PD_Peptides_95pct',
           'Gene_Length_bp', 'GC_Content_Pct', 'Mean_eQTL_SNP_Count']
t1 = covar[t1_cols].copy()

# Source ordering: Candidates first (Pull-down then Literature), then Non-Candidates, then T2DM
group_order = {'Candidate': 0, 'NonCandidate': 1, 'T2DM_Control': 2}
t1['_order'] = t1['Group'].map(group_order)
t1 = t1.sort_values(['_order', 'Source', 'Gene']).drop(columns='_order')

t1.to_csv(os.path.join(OUT, 'TableS1_gene_list.csv'), index=False)
print(f"  → Saved: TableS1_gene_list.csv ({len(t1)} genes)")
print(f"    Groups: {t1['Group'].value_counts().to_dict()}")
print(f"    Sources: {t1['Source'].value_counts().to_dict()}")

# ─────────────────────────────────────────────
# Table S2: GTEx Baseline TWAS Results
# ─────────────────────────────────────────────
print("\nGenerating Table S2: GTEx Baseline TWAS Results...")

# Read all GTEx tissue files
tissues_data = []
for tissue_label, prefix in [('Nerve_Tibial', 'Nerve_Tibial'), ('Whole_Blood', 'Whole_Blood')]:
    for trait in ['DR', 'DN', 'DPN']:
        fn = f'gtex_{prefix}_{trait}.csv'
        df = pd.read_csv(os.path.join(DATA, fn))
        df = df.rename(columns={
            'gene': 'Gene',
            'zscore': f'Z_{tissue_label}_{trait}',
            'pvalue': f'P_{tissue_label}_{trait}',
            'n_snps_model': f'SNPs_Model_{tissue_label}_{trait}',
        })
        df = df[['Gene', f'Z_{tissue_label}_{trait}', f'P_{tissue_label}_{trait}', f'SNPs_Model_{tissue_label}_{trait}']]
        tissues_data.append(df)

# Merge all tissue-trait columns
t2 = tissues_data[0]
for df in tissues_data[1:]:
    t2 = t2.merge(df, on='Gene', how='outer')

# Add multi-tissue integrated results
stouffer = pd.read_csv(os.path.join(DATA, 'gtex_stouffer_integrated.csv'))
for trait in ['DR', 'DN', 'DPN']:
    sub = stouffer[stouffer['Trait'] == trait][['Gene', 'Z_Multi', 'P_Multi', 'FDR_q']].copy()
    sub = sub.rename(columns={'Z_Multi': f'Z_Multi_{trait}', 'P_Multi': f'P_Multi_{trait}', 'FDR_q': f'FDR_Multi_{trait}'})
    t2 = t2.merge(sub, on='Gene', how='outer')

# Add ACAT-O results
acat = pd.read_csv(os.path.join(DATA, 'gtex_acat_o_results.csv'))
for trait in ['DR', 'DN', 'DPN']:
    sub = acat[acat['Trait'] == trait][['Gene', 'P_ACAT_O', 'FDR_q']].copy()
    sub = sub.rename(columns={'P_ACAT_O': f'P_ACAT_O_{trait}', 'FDR_q': f'FDR_ACAT_O_{trait}'})
    t2 = t2.merge(sub, on='Gene', how='outer')

# Add group info
group_map = dict(zip(covar['Gene'], covar['Group']))
t2['Group'] = t2['Gene'].map(group_map)
t2 = t2[['Gene', 'Group'] + [c for c in t2.columns if c not in ['Gene', 'Group']]]

# Sort: candidates first
group_order_s2 = {'Candidate': 0, 'NonCandidate': 1, 'T2DM_Control': 2}
t2['_order'] = t2['Group'].map(group_order_s2)
t2 = t2.sort_values(['_order', 'Gene']).drop(columns='_order')

t2.to_csv(os.path.join(OUT, 'TableS2_gtex_baseline_results.csv'), index=False)
print(f"  → Saved: TableS2_gtex_baseline_results.csv ({len(t2)} genes, {len(t2.columns)} columns)")
print(f"    Columns: {list(t2.columns)}")

print("\nDone!")
