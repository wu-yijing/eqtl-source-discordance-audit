# Processed data — eQTL weight-source audit of TWAS gene candidacy

Supporting processed data for:

> **"eQTL weight-source choice reshapes TWAS gene candidacy: a two-axis dual-source audit with disease-agnostic calibration and a genome-wide benchmark"** (BMC Genomics submission)

> **Scope of this file.** It describes what is in `data/`. It is **not** the authoritative statement of how the archive maps onto the manuscript — that is `../ARCHIVE_NOTE.md`, together with `processed_officialZ/_PROVENANCE.json`. The inventory of `processed/` further down is retained **as a historical record** of the early-generation (pre-correction) intermediate layer.

## Which layer is authoritative

| Layer | Status | Use it for |
|---|---|---|
| `processed_officialZ/` | ✅ **Authoritative** — official MetaXcan v0.8.1 recompute | Every headline value in the manuscript: GTEx and eQTLGen Z-scores, FDR calls, denominators, cross-cohort values and the genome-wide SCZ arms |
| `processed/` | ⚠️ **Superseded** — pre-correction (missing S-PrediXcan σᵢ expression-variance factor and a PLINK 2-bit decoding defect) | Nothing in the manuscript, except the "archived working table" that Additional file 1: Table S10 names for the four anchor-set proteins |

The deprecation record is `processed/_PROVENANCE.json`; a narrative version is `processed/_DEPRECATED_勿用_修正前数据_20260917.md`. The two defects are disclosed in the manuscript's Methods, and the earlier routine is retained only as an equivalence cross-check (residual |ΔZ| ≤ 3 × 10⁻⁸).

## The 104-gene testbed

- **Groups:** 30 HOTAIR-interactome candidate + 44 non-candidate + 30 T2DM control genes
- **Phenotypes:** diabetic retinopathy (DR), diabetic nephropathy (DN), diabetic peripheral neuropathy (DPN)
- **eQTL sources:** GTEx v8 MASHR (Nerve_Tibial, Whole_Blood; multi-tissue ACAT-O and Stouffer integration) and eQTLGen phase I whole blood (N = 31,684)
- **GWAS input:** FinnGen R13 (DR/DN/DPN); UK Biobank GCST90043640 for the cross-cohort check; PGC3 wave-3 schizophrenia for the genome-wide benchmark

## Authoritative files (`processed_officialZ/`)

| File | Rows | Content |
|---|---|---|
| `gtex_official_Z.csv` | 222 | GTEx v8 MASHR multi-tissue per-gene Z, P and FDR q for 74 testbed genes × 3 phenotypes (= Additional file 1: Table S2) |
| `eqtlgen_official_Z.csv` | 288 | eQTLGen whole-blood per-gene Z, P, BH q and model-SNP counts for 96 genes (27 candidate + 25 non-candidate + 17 T2DM + 27 housekeeping); the 69-gene universe of Additional file 1: Table S17 is a subset |
| `primary_arm_96pairs_official.csv` | 96 | The pre-specified primary comparison, both source Z-scores per gene–phenotype pair (= Additional file 1: Table S12) |
| `crosscohort_TableS4_official.csv` | 4 | Cross-cohort values (= Additional file 1: Table S4a) |
| `gene_groups_TableS1_official.csv` | 104 | Gene list with group assignments (= Additional file 1: Table S1) |
| `scz_z_4arm_official.csv` | 15,875 | Genome-wide PGC3 SCZ four-arm Z-scores (8,315-gene complete-case set) |
| `_PROVENANCE.json` | — | Which layer is authoritative, and why |

## Historical inventory of `processed/` (pre-correction; retained for provenance)

| File | Rows | Description as originally generated |
|---|---|---|
| `eqtlgen_spredixcan_results.csv` | 219 | eQTLGen S-PrediXcan full results |
| `eqtlgen_spredixcan_harmonized_results.csv` | 309 | Three-way allele-harmonised eQTLGen results |
| `eqtlgen_vs_gtex_comparison.csv` | 102 | Paired GTEx vs eQTLGen Z-scores with direction indicators |
| `covariate_matrix.csv` | 104 | Gene-level covariates (group, protein source, gene length, GC%, eQTL SNP counts) |
| `mahalanobis_matched_pairs.csv` | 60 | 30 Mahalanobis-matched candidate–control pairs with covariate values |
| `gtex_acat_o_results.csv` | 237 | GTEx multi-tissue ACAT-O results for 79 genes — **the "archived working table" named in Additional file 1: Table S10** |
| `gtex_stouffer_integrated.csv` | 237 | Same gene universe, Stouffer integration |
| `enrichment_comparison.csv`, `enrichment_comparison_harmonized.csv` | — | FDR enrichment rates by group, phenotype and source |
| `candidate_comparison_DR.csv`, `viz_z_distribution.csv`, `layer_analysis.csv` | — | Visualisation and stratification inputs |
| `gtex_Nerve_Tibial_*.csv`, `gtex_Whole_Blood_*.csv` | 62–67 each | Single-tissue per-gene results (early generation) |

> **Note on the matched pairs.** `mahalanobis_matched_pairs.csv` is listed among the superseded files because it was written in the same pre-correction run, but its contents are **covariate-only** (gene length, GC content, eQTL SNP counts, matched subclass). The two S-PrediXcan defects cannot propagate to covariates: the 30 matched pairs and every numeric covariate value are **identical** to those in Additional file 1: Table S3. The matched enrichment contrasts of Additional file 1: Table S25 are reproduced by combining this file with the *authoritative* Z tables `processed_officialZ/gtex_official_Z.csv` and `processed_officialZ/eqtlgen_official_Z.csv`; see `../ARCHIVE_NOTE.md`.

## Data sources (all publicly available)

| Dataset | Access |
|---|---|
| FinnGen R13 (2024 release) | https://www.finngen.fi/en/access_results |
| GTEx v8 MASHR weights (models of 21 Jul 2021) | https://predictdb.org/ |
| eQTLGen phase I cis-eQTL statistics | https://www.eqtlgen.org/ |
| 1000 Genomes phase 3 EUR LD panel | https://www.internationalgenome.org/ |
| UK Biobank DR, GWAS Catalog GCST90043640 | https://www.ebi.ac.uk/gwas/studies/GCST90043640 |
| DN atlas, Sakaue et al. 2021 (ebi-a-GCST90018832) | https://www.ebi.ac.uk/gwas/studies/GCST90018832 |
| PGC3 schizophrenia wave 3 | figshare, DOI 10.6084/m9.figshare.19426775 |
| HRT Atlas v1.0 housekeeping universe | https://www.hrtatlas.com/ |

## License

The code **and** the processed data in this repository are released under the **MIT licence** (see `../LICENSE`), matching the manuscript's Data availability statement.

*(Corrected in v1.0.1: earlier revisions of this file carried a superseded manuscript title, a GigaDB dataset framing, a CC0-1.0 licence statement and several pre-correction row counts. See `../ARCHIVE_NOTE.md`.)*
