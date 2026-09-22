# GigaDB Dataset: TWAS eQTL Source Confounding Evaluation

## Supporting data for: "eQTL source selection systematically biases TWAS cross-population replication"

### Dataset Overview

This dataset contains processed S-PrediXcan TWAS results comparing two eQTL weight sources (GTEx v8 MASHR multi-tissue vs. eQTLGen whole blood) for 104 genes across three diabetic microvascular complications, using FinnGen R13 GWAS data.

**Genes:** 30 HOTAIR candidate + 44 non-candidate + 30 T2DM controls  
**Phenotypes:** DR (diabetic retinopathy), DN (diabetic nephropathy), DPN (diabetic peripheral neuropathy)  
**eQTL sources:** GTEx v8 (Nerve_Tibial + Whole_Blood, N≈600) and eQTLGen (Whole Blood, N=31,684)  
**GWAS source:** FinnGen Release 13

### File Inventory

| File | Rows | Columns | Size | Description |
|------|------|---------|------|-------------|
| `eqtlgen_spredixcan_results.csv` | 219 | 8 | 18 KB | eQTLGen S-PrediXcan full results: Z-scores and P-values for all gene-phenotype pairs |
| `eqtlgen_vs_gtex_comparison.csv` | 102 | 6 | 5 KB | Paired GTEx vs eQTLGen Z-scores with direction consistency indicators |
| `covariate_matrix.csv` | 114 | 14 | 9 KB | Gene-level covariates: group assignment, protein source, length, GC%, eQTL SNP counts |
| `mahalanobis_matched_pairs.csv` | 60 | 8 | 3.6 KB | 30 Mahalanobis-matched candidate-control gene pairs with covariate values (updated: 2026-06-30, 44 high-confidence non-candidate pool) |
| `enrichment_comparison.csv` | 6 | 10 | 618 B | FDR enrichment rates by group, phenotype, and eQTL source |
| `candidate_comparison_DR.csv` | 28 | 9 | 2.4 KB | Candidate gene eQTLGen Z-scores for DR phenotype with GTEx cross-reference |
| `viz_z_distribution.csv` | 126 | 7 | 7.5 KB | Z-score distribution data for density and scatter visualization |
| `layer_analysis.csv` | 6 | 9 | 319 B | Pull-down vs literature source stratification results |

### Column Descriptions

#### eqtlgen_spredixcan_results.csv
- `Gene` - Gene symbol
- `Trait` - Diabetic complication (DR/DN/DPN)
- `Tissue` - eQTL source (eQTLGen_Whole_Blood)
- `Group` - Gene group (Candidate/NonCandidate/T2DM_Control)
- `Z_eQTLGen` - TWAS Z-score using eQTLGen weights
- `P_eQTLGen` - TWAS P-value
- `n_SNPs_Model` - Number of SNPs in the predictive model
- `n_SNPs_Matched` - Number of SNPs matched to GWAS summary statistics

#### eqtlgen_vs_gtex_comparison.csv
- `Gene` - Gene symbol
- `Trait` - Diabetic complication (DR/DN/DPN)
- `Group` - Gene group
- `Z_eQTLGen` - TWAS Z-score using eQTLGen weights
- `Z_GTEx` - TWAS Z-score using GTEx v8 multi-tissue weights
- `Same_Direction` - Boolean: whether both Z-scores have the same sign

#### covariate_matrix.csv
- `Gene` - Gene symbol
- `Group` - Gene group (Candidate/NonCandidate/T2DM_Control)
- `Source` - Protein identification method (Pull-down/Literature_curation)
- `PullDown_Unused` - Protein confidence score (Unused)
- `PullDown_Cov_pct` - Protein sequence coverage percentage
- `PullDown_Peptides95` - Number of peptides at 95% confidence
- `Length_bp` - Gene length in base pairs
- `GC_pct` - Gene GC content percentage
- `Max_eQTL_SNPs` / `Min_eQTL_SNPs` / `eQTL_SNPs_Max` / `eQTL_SNPs_Min` / `eQTL_SNPs_Mean` / `eQTL_SNPs_Mean_Num` - eQTL SNP count statistics

### Data Sources (publicly available)

| Dataset | Access |
|---------|--------|
| FinnGen R13 | https://www.finngen.fi/ |
| GTEx v8 MASHR weights | https://gtexportal.org/ |
| eQTLGen | https://www.eqtlgen.org/ |
| 1000 Genomes EUR LD | https://www.internationalgenome.org/ |
| UK Biobank DR (GCST90043640) | https://gwas.mrcieu.ac.uk/ |
| Xue et al. 2022 DR | ⚠️ **Not resolvable** — the cited ID (`ieu-b-4803`) is absent from IEU OpenGWAS (verified 2026-09-16); the publication's real DR subset (GCST90134546, 1,652 cases / 60,577 controls) has no public summary statistics on the GWAS Catalog FTP. See `analysis_reports/数据集编号溯源_Xue2022_ieu-b-4803.md`. |
| Sakaue et al. 2021 DN | https://gwas.mrcieu.ac.uk/ |

### Analysis Code

All analysis scripts are available at:
https://github.com/wu-yijing/eqtl-source-discordance-audit
(MIT License)

### License

This dataset is licensed under CC0 1.0 Universal.
