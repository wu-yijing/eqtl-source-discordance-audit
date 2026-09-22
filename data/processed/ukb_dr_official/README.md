# Official MetaXcan results for the cross-cohort DR arm

All Z-scores were produced with the unmodified official `SPrediXcan.py` from MetaXcan v0.8.1
(Python 3.12 + numpy 1.26; `PYTHONWARNINGS=ignore` is set in the caller only, the source is untouched),
on the 1000 Genomes EUR panel covariance built in this project, with the GWAS supplied using its own
effect/non-effect allele columns so that MetaXcan performs the allele alignment itself
(`metax/misc/GWASAndModels.py:align_data_to_alleles`).

| File | GWAS | Weight source |
|---|---|---|
| `finngen_r13_dr_eqtlgen_official.csv` | FinnGen R13 DR | eQTLGen whole blood |
| `ukb_gcst90043640_eqtlgen_official.csv` | UKB GCST90043640 | eQTLGen whole blood |
| `ukb_gcst90043640_gtex_Nerve_Tibial_official.csv.gz` | UKB GCST90043640 | GTEx v8 MASHR Nerve_Tibial |
| `ukb_gcst90043640_gtex_Whole_Blood_official.csv.gz` | UKB GCST90043640 | GTEx v8 MASHR Whole_Blood |

`RNH1_official_metaxcan_Z.csv` extracts the RNH1 row of each run.

Protocol note: the FinnGen × eQTLGen run was repeated from the raw FinnGen allele columns
(`finngen_r13_dr_eqtlgen_official.csv`, RNH1 Z = +2.3091, 726/898 model SNPs) to validate that the
official binary performs allele alignment; the pre-aligned input used for the arm-level recompute
returns Z = +2.3064 (728/898). The two differ by 0.003, attributable to duplicate-rsID handling.

Retraction note: the formerly used second cohort `ieu-b-4803` is **not resolvable** in IEU OpenGWAS
(verified 2026-09-16 against the full 50,057-dataset catalogue) and is no longer part of any analysis.
See `analysis_reports/数据集编号溯源_Xue2022_ieu-b-4803.md`.
