# eQTL weight-source choice reshapes TWAS gene candidacy — a two-axis dual-source audit

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Repository scope

This repository holds **analysis code and processed data only** for:

> **"eQTL weight-source choice reshapes TWAS gene candidacy: a two-axis dual-source audit with disease-agnostic calibration and a genome-wide benchmark"**

Status: code & data release accompanying the manuscript (BMC Genomics submission).

**Not in this repository — by design:**

| Excluded | Where it lives |
|---|---|
| All figures and figure captions | The manuscript and *Additional file 1* (journal supplementary files) |
| All tables and table legends | The manuscript and *Additional file 1* |
| The manuscript, cover letter and supplementary `.docx` files | Journal submission only |
| Raw RNA pull-down / LC–MS/MS spectra and search results | ProteomeXchange Consortium via iProX (identifier quoted in the manuscript) |
| Internal working notes and drafting reports | Not distributed (`analysis_reports/` is `.gitignore`d) |

`figures/` and `tables/` are listed in `.gitignore` and are **not** tracked. `figure_scripts_officialZ_20260917/`
writes figure output into `<repo>/figures/` at run time — that directory is a runtime output, not repository content.

### Core Question

By how much does the choice of eQTL weight source (GTEx v8 tissue-specific vs. eQTLGen large-scale whole blood) perturb the conclusions a TWAS actually reports — Z-score sign, candidate gene rank and FDR gene candidacy — and how much of that perturbation is attributable to the resource/sample-size axis rather than the tissue-context axis?

> **Terminology.** Throughout this repository and the current manuscript, *eQTL-source discordance* denotes the systematic perturbation of TWAS gene-candidacy conclusions induced by the weight source when GWAS input and analytic pipeline are held fixed. It is operational and **not** a claim of confounding in the causal sense.

## Key Findings

> **Provenance.** The values below are the **current manuscript values** (BMC Genomics submission), recomputed with the **unmodified official MetaXcan v0.8.1 binary**. `data/processed_officialZ/_PROVENANCE.json` records the archived inputs, and `figure_scripts_officialZ_20260917/` is the pipeline that reproduces the manuscript's figure values from them.

> **Two corrections these numbers rest on.** (i) **Allele harmonization**: all eQTLGen-arm values were recomputed after three-way allele harmonization (eQTLGen assessed/other alleles ↔ FinnGen R13 ref/alt ↔ 1000G EUR A2) against the official eQTLGen FDR < 0.05 weight catalogue; the earlier archived eQTLGen export carried an allele-alignment defect (24.1% of Z-scores flipped sign). (ii) **Estimator correction**: the earlier in-house S-PrediXcan implementation omitted the σᵢ expression-variance factor and mis-decoded the PLINK 2-bit genotype encoding. Both were fixed and every Z-score recomputed with the official MetaXcan v0.8.1 binary; this is what removed the previously reported extreme values (RNH1 DR Z fell from 13.32 to 2.31, TUBB from 48.52 to 11.89, CKAP4 from 4.75 to 0.99).

> **Definitions and denominators.** FDR enrichment rates are Benjamini–Hochberg-corrected **within each group × phenotype × weight-source stratum** and are the proportion of testable genes with q < 0.05, compared between groups by Fisher's exact test. Two eQTLGen gene universes must not be conflated: the **61-gene harmonised arm** (24 candidate, 24 non-candidate, 13 T2DM control), which carries the model-SNP statistics, and the **69-gene denominator universe** (27 candidate, 25 non-candidate, 17 T2DM control; 207 gene–phenotype pairs), which supplies the eQTLGen rates of Table 2 and the per-gene FDR calls of Additional file 1: Tables S17 and S18. Of the 104 testbed genes, 37 have no eQTLGen model in the official FDR < 0.05 catalogue and 67 do, of which 61 return a valid S-PrediXcan statistic (Additional file 1: Table S10). Gene-level Z-scores and the arm denominators are in `data/processed_officialZ/`.

| Metric | Value | Notes |
|--------|-------|-------|
| Total genes analysed (testbed) | 104 = 30 candidate + 44 non-candidate + 30 T2DM control | One row per gene in the covariate matrix (Additional file 1: Table S1). The three disease-agnostic control **layers** (30 genes each; Table 1a) are a separate set, not a fourth testbed group. |
| Spearman ρ, GTEx Whole_Blood vs eQTLGen (primary arm) | 0.39 (Fisher-z 95% CI 0.21–0.55; naive P = 8.8 × 10⁻⁵) | 96 gene–phenotype pairs from 32 genes. The naive P is anticonservative because pairs cluster within genes; cluster-aware alternatives give SE = 0.125 (two-sided P = 0.004), jackknife SE = 0.137 (P = 0.008) and permutation P < 0.001 (Additional file 1: Table S16). |
| Direction consistency (primary arm) | 68.8% (66/96); gene-cluster bootstrap 95% CI 58.3–79.2% | DR 71.9% (23/32), DN 65.6% (21/32), DPN 68.8% (22/32) (Table 3a). Not separable from the within-arm empirical null (63.6%, 35/55) and 27.9 percentage points below the 96.7% within-source reference (95% CI 17.5–38.4) — a monotone transform of ρ, so the shortfall is reported as a bound, not an effect size. |
| Direction consistency by gene-subset definition (Table 3c) | 67.6% (69/102) / 68.8% (66/96) / 71.1% (113/159) | Both-source anchor set (34 genes, ρ = 0.373) / primary arm (32 genes, ρ = 0.390) / panel-only all-testbed (53 genes, ρ = 0.490). Across the three whole-testbed definitions the rate spans 67.6–71.1% and ρ spans 0.373–0.490; the candidate-only definitions reach 74.6% (47/63, panel-only) and 74.4% (58/78, dual). |
| GTEx FDR enrichment (candidate arm, all phenotypes) | 2.4% (2/84) | GTEx v8 MASHR with multi-tissue ACAT-O integration. Same-source comparators: non-candidates 1.2% (1/81), T2DM controls 1.8% (1/57) — candidates are not enriched above either (Fisher's exact P = 0.24 and 0.61). Per-phenotype candidate rates are 0.0–3.6% (Tables 1b and 2). |
| eQTLGen FDR enrichment (candidate arm, all phenotypes) | 6.2% (5/81) | Non-candidates 0.0% (0/75) and T2DM controls 0.0% (0/51); no significant excess (Fisher's exact P = 0.21 and 0.49). After Mahalanobis matching, 6.2% (5/81) against 0.0% (0/57) in matched controls (P = 0.077); the GTEx arm is indistinguishable (2.4% vs 1.7%; P = 1.00). No group exceeded 7.4%, below the ≥8.0-pp (GTEx) and ≥13.0-pp (eQTLGen) 80%-power boundaries. |
| RNH1 GTEx Z-score (DR, Nerve_Tibial) | +2.67 | S-PrediXcan; DN +1.70, DPN +1.63; ACAT-O multi-tissue q = 0.28, 0.67, 0.62 (none FDR-significant). The only candidate pair reaching \|Z\| ≥ 1.96 in both sources (Additional file 1: Table S17). |
| RNH1 eQTLGen Z-score (DR) | +2.31 | N = 31,684; P = 0.021. Same sign as the GTEx estimate — the source-stable index gene. |
| Cross-cohort heterogeneity (RNH1, DR; FinnGen R13 vs UK Biobank GCST90043640, both eQTLGen weights) | Z scale (primary): pooled +1.51 (SE 0.79; P = 0.056), Q = 1.26, I² = 20.6%, 95% prediction interval −0.33 to +3.36 (spans zero). √N_e-weighted re-merge: +2.39 (P = 0.017), Q = 0.13, I² = 0% | FinnGen Z = +2.31 (N_e = 49,304), UKB Z = +0.72 (N_e = 1,231; 308 cases), both direction-consistent. With k = 2 and an ≈40-fold N_e gap the heterogeneity statistics are convention-dependent (the same weights on the raw Z scale give +2.09, Q = 76.6, I² = 98.7%), so the between-cohort contrast, not Q or I², carries the evidence (Additional file 1: Table S4). |
| Genome-wide SCZ replication (PGC3 wave3 European) | ρ = 0.469 panel-only / 0.420 tissue-only / 0.447 dual (n = 8,315) | Complete-case **8,315** genes (80.3%) of a **10,357**-gene pool; 9,048 (87.4%) carry a GTEx Whole_Blood model and 8,890 a Z in both GTEx tissues. The plotted tissue-axis estimate uses n = 8,890, the four-arm decomposition n = 8,315, and both give ρ = 0.420 to three decimals. Arm-wise direction consistency 67.2 / 66.8 / 66.2% (pairwise gaps ≤ 1.0 pp); dual-arm discordance 33.8% (Additional file 1: Table S15). |

## Genome-Wide SCZ 2×2 Decomposition

A dual-source 2×2 decomposition was applied to schizophrenia (SCZ, PGC3 wave3 European, 53,386 cases / 77,258 controls) to test whether the source/panel vs tissue-context structure generalizes beyond the metabolic discovery cohort.

- **Method**: the SCZ statistics were recomputed with the **official MetaXcan v0.8.1 binary**, holding the eQTL weights, tissues, SNP filters and LD reference identical to the discovery analysis. The four-arm Z matrix is archived as `data/processed_officialZ/scz_z_4arm_official.csv`.
- **Denominators (two are in use)**: a 10,357-gene pool with a valid model, of which 9,048 (87.4%) carry a GTEx Whole_Blood model and 8,890 a Z in **both** GTEx tissues. The four-arm decomposition requires a finite Z in all four arms and therefore rests on the **8,315-gene complete-case set (80.3% of the pool)**. The plotted tissue-axis estimate uses n = 8,890 and the decomposition n = 8,315; the two subsets agree to three decimals (ρ = 0.420).
- **Results (complete-case n = 8,315)**:
  - Panel-only analogue (eQTLGen vs GTEx Whole_Blood): ρ = **+0.469**, 67.2% direction-consistent
  - Tissue-only analogue (GTEx Whole_Blood vs GTEx Nerve_Tibial): ρ = **+0.420**, 66.8%
  - Dual analogue (eQTLGen vs GTEx multi-tissue): ρ = **+0.447**, 66.2%, reversal rate 33.8% (against 32.8% in the testbed dual arm)
  - The three arms are indistinguishable on direction consistency (every pairwise gap ≤ 1.0 percentage point, against a binomial interval width of ≈±1 point), while on ρ the dual arm lies between the two single-axis arms — the ρ-based axis ordering is **trait-dependent**, not a fixed "source always dominates" rule.
  - Formal axis contrast: Δρ(dual − tissue) = +0.0265 (95% CI +0.0002 to +0.0524; bootstrap P = 0.048) and Δρ(dual − panel) = −0.0226 (95% CI −0.0355 to −0.0094; P = 0.0006). Measured against the post hoc 0.10 margin the axes are reported as **indistinguishable at this sample size**, not as equivalent (the margin is 3.8× the observed difference).
- **Sparsity**: GTEx v8 MASHR models are sparse (median 2 model SNPs per gene, IQR 1–2) whereas the eQTLGen cis panel is dense, and the tissue-axis contrast is not framework-dependent (ρ = 0.525 under elastic net vs 0.499 under MASHR on the 4,098-gene four-way universe). Stratifying the SCZ contrast by GTEx model size does not establish sparsity as a sufficient explanation of the tissue axis (tissue-only minus dual, +0.54 percentage points, 95% CI −0.79 to +1.86, P = 0.43; Additional file 1: Table S15).
- **Authoritative values**: `data/processed_officialZ/scz_z_4arm_official.csv` together with Additional file 1 (Tables S15 and S16).

## Quick Start (Docker — Recommended)

```bash
# 1. Build the Docker image
docker build -t eqtl-discordance-audit .

# 2. Run the full reproducibility pipeline
docker run --rm -v $(pwd)/output:/app/output eqtl-discordance-audit

# 3. View outputs
open output/figs/    # figures regenerated by the container pipeline (PDF + PNG)
open output/logs/    # analysis provenance logs
```

> ⚠️ **What the container actually runs.** `run_all.sh` executes the **early-generation** downstream scripts in `scripts/python/` (they write `figs/`, which is then copied to `output/`). The container therefore reproduces the environment and that earlier pipeline — **not** the manuscript's current figure values. To regenerate the current figures use `figure_scripts_officialZ_20260917/`; its step 00 needs *Additional file 1* (the journal supplementary file, not distributed here), so it cannot run in the container without you supplying that file.

### What the Docker image contains
- **Python 3.13 + R 4.5.2** with all dependencies (pandas, scipy, statsmodels, MatchIt, etc.)
- **MetaXcan** (for S-PrediXcan, run separately with `run_spredixcan.sh`)
- **Downstream analysis scripts**
- **Processed data** needed to reproduce results from intermediate results

## Local Setup (Conda)

```bash
conda env create -f environment.yml
conda activate twas-eqtl

# Run downstream analysis
python scripts/python/04_generate_all_figures.py
```

## Reproducing the Full Pipeline (Optional)

To reproduce from raw GWAS data (requires large input files):

```bash
docker run --rm \
  -v /path/to/finngen_data:/input/finngen \
  -v /path/to/gtex_models:/input/gtex \
  -v /path/to/eqtlgen_weights:/input/eqtlgen \
  -v /path/to/1000g_ref:/input/1000g \
  -v $(pwd)/output:/app/output \
  eqtl-discordance-audit \
  bash /app/run_spredixcan.sh
```

## Figure pipeline (`figure_scripts_officialZ_20260917/`)

This directory is the **only** pipeline that reproduces the manuscript's current figure values: it reads
`data/processed_officialZ/` (official MetaXcan v0.8.1 Z recompute) and *Additional file 1*'s tables.

All paths go through `paths_config.py` — run `python paths_config.py` to print the resolved input/output
paths, and override with `TWAS_REPO` / `TWAS_DATA_Z` / `FIG_OUT_MAIN` / `FIG_OUT_SUPP` / `AF1_DOCX` / `FIG_RESULTS`.
**Additional file 1 is not distributed with this repository** (it is the journal's supplementary file);
download it and point `AF1_DOCX` at it. See that directory's `README.md`.

Figure **output** goes to `<TWAS_REPO>/figures/`, which is `.gitignore`d — the figure files themselves are
not part of this repository.

### Earlier script generations (not distributed)

| Generation | Where it was | Why it is not here |
|---|---|---|
| Five-script figure set (pre-2026-09-17) | `figure_scripts/` | Read the **pre-correction** data layer `data/processed/` (missing the σᵢ expression-variance factor, with a PLINK 2-bit decoding defect), so it reproduces *retired* values — e.g. RNH1 DR Z = 13.318 vs 2.3064, TUBB = 48.516 vs 11.8932, CKAP4 = 4.7549 vs 0.9874 (×4–×6). |
| Early supplementary-figure scheme | `FigureS1–FigureS8` | Does not correspond to the current *Additional file 1* figure numbering. |
| In-house SCZ TWAS estimator | `scz_replication/` | Implements `Z = (w′R⁻¹z)/√(w′R⁻¹w)` rather than the official MetaXcan v0.8.1 binary. Its snapshot (n = 2,511) is not the manuscript's value. Its cluster-robustness re-analysis is deliberately still shipped, under `scripts/python/s1_cluster_robustness/`, because it documents that earlier run. |

## Repository Structure

```
├── README.md                       # this file
├── LICENSE (MIT)  .gitignore  .dockerignore  .zenodo.json
├── Dockerfile  environment.yml  requirements.txt  renv.lock
├── run_all.sh                      # pipeline entry point (Docker ENTRYPOINT)
├── run_spredixcan.sh               # S-PrediXcan reference (needs external data)
├── run_three_figs.py
│
├── figure_scripts_officialZ_20260917/   # CURRENT figure pipeline (writes to a gitignored figures/)
│   ├── paths_config.py             #   single path entry point (env-overridable)
│   ├── 00…08, 10, 11 *.py          #   00 = data layer; 01/02/04/06/08/10 = figures; 05/07/11 = checks
│   ├── m15_positive_control.json   #   input behind Fig. S6
│   └── fig4_bootstrap_officialZ.json   # output of 06 (B = 5,000, seed 20260915)
│
├── data/
│   ├── processed_officialZ/        # AUTHORITATIVE Z / denominator layer (built from Additional file 1)
│   │   └── _PROVENANCE.json
│   ├── processed/                  # early-generation intermediates — DEPRECATED for Z-based results
│   │                               #   (file-level notices inside; see data/README.md)
│   ├── hk_reselect_20260830/       # control-layer (housekeeping) selection data
│   └── raw/                        # NOT tracked — raw public GWAS/eQTL data, too large for git
│
├── scripts/
│   ├── python/                     # EARLY-generation pipeline (provenance only)
│   │                               #   reads data/processed/, writes a runtime figs/; not the
│   │                               #   current figure pipeline — see its README_DEPRECATED.md
│   └── R/run_mahalanobis_matching.R
├── analyses/                       # control-layer scripts + logs/ (run provenance)
├── audit_notes/                    # dated audit notes shipped with the repo — see audit_notes/README.md
│                                   #   NB: analysis_reports/ is a local-only working directory, not tracked
```

**Which part is current.** The Z-scores, denominators and cross-cohort values are authoritative in
`data/processed_officialZ/`, and the only pipeline that reproduces the manuscript's figure values is
`figure_scripts_officialZ_20260917/` (start with `python paths_config.py`). Older generations are kept for
provenance and marked either in-file, in `scripts/python/` (early-generation pipeline) or in
`data/processed/` (pre-correction intermediates) — none of them should be used to reproduce the paper.
A directory named `figs/` is **not** part of the repository: `run_all.sh` creates it at run time as its
figure output, then copies it (with logs and processed data) into `output/`.

## Cluster-Robustness Re-analysis of Direction Consistency (S1)

> ⚠️ **Status.** This section describes the corrected Z layer. `scripts/python/s1_cluster_robustness/`
> still documents the earlier generation and should be read as provenance for it, not as the current
> analysis. The values below are the manuscript's — Additional file 1: Table S16 for the cluster-aware
> procedures and Table 3a for the arm rates.

The primary between-source comparison comprises 96 gene–phenotype pairs that derive from only **32 genes**
(three phenotypes each), so the pairs are not independent observations. All four cluster-aware procedures
reported by the manuscript are computed on the rank correlation (ρ = 0.39 on the 96 pairs):

| Procedure | Value |
|---|---|
| Naive t test on the rank correlation | t = 4.10 (df 94), P = 8.8 × 10⁻⁵ |
| Cluster-robust sandwich standard error | SE = 0.125 — one-sided P = 0.002 (two-sided P = 0.004, df 31) |
| Delete-one-gene jackknife standard error | SE = 0.137 — one-sided P = 0.004 (two-sided P = 0.008, df 31) |
| Gene-label permutation (B = 10,000) | null 2.5–97.5 percentiles −0.22 to 0.23; P < 0.001 |
| Gene-level cluster bootstrap (B = 10,000) | ρ 95% CI 0.12–0.62; direction consistency 95% CI 58.3–79.2% |

Resampling detail: all four resampling procedures use B = 10,000 draws with **seed 20260915**
(Additional file 1: Table S16, note).

All four procedures keep the association significant, so the independent-pair naive test is not materially
anti-conservative. The manuscript nevertheless reports the **gene-cluster bootstrap interval as the primary
interval** for the direction-consistency rate (58.3–79.2%) rather than the independent-pair Clopper–Pearson
interval (58.5–77.8%; Additional file 1: Table S20, min|Z| = 0), and reports the rate descriptively, because
it is not separable from the within-arm empirical null (63.6%, 35/55).

**Repo-level diagnostic — not reported in the manuscript.** Recomputing the one-way ICC of the per-gene
concordance proportions on the current primary arm
(`data/processed_officialZ/primary_arm_96pairs_official.csv`; 32 genes × 3 pairs) gives **ICC = 0.000**
(the raw estimate is −0.004 before clipping at zero) and **design effect = 1.000** — between-gene variance
does not exceed within-gene variance. Clustering of the binary concordance outcome within genes is therefore
negligible. This is offered for reproducibility and is not part of the manuscript's evidence.

**Genome-wide SCZ arm.** With a single trait each gene contributes **exactly one comparison** (8,315
complete-case genes drawn from a 10,357-gene pool), so the design effect is 1 by construction and no cluster
correction is required. The arm correlations are ρ = +0.469 (panel-only), +0.420 (tissue-only) and +0.447
(dual), with direction consistency 67.2 / 66.8 / 66.2%.

## Data Sources

All GWAS and eQTL summary statistics used in this study are from **publicly available sources**:

| Dataset | Source | Access |
|---------|--------|--------|
| FinnGen R13 (DR/DN/DPN) | [FinnGen](https://www.finngen.fi/) | Release 13 |
| GTEx v8 MASHR weights | [GTEx Portal](https://gtexportal.org/) | Public access |
| eQTLGen whole blood cis-eQTL | [eQTLGen](https://www.eqtlgen.org/) | Public access |
| 1000 Genomes EUR LD | [1000 Genomes](https://www.internationalgenome.org/) | Public access |
| UK Biobank DR (GCST90043640) | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/) | Public access |
| ~~UK Biobank DR (Xue et al. 2022; ieu-b-4803, 14,147 cases/322,390 controls)~~ | ⚠️ **Withdrawn 2026-09-16** — this accession is not present in IEU OpenGWAS (verified against the full 50,057-dataset catalogue; the neighbouring IDs `ieu-b-4795`–`ieu-b-4807` are likewise absent), and no local or repository copy exists; the Z value formerly attributed to it had no traceable data or code. The real DR subset of that publication (PMID 35841873) is GCST90134546 (1,652 cases / 60,577 controls) and has no public summary statistics. **Replaced by the already-used UKB DR dataset GCST90043640 (308 cases / 456,040 controls, re-analysed with eQTLGen weights).** |
| Cai et al. 2026 UKB T2D-DR | (in preparation) | Apply to UKB |
| Sakaue et al. 2021 DN (ebi-a-GCST90018832) | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/) | Public access |

## Requirements

### Python
```bash
pip install -r requirements.txt
```

Or with conda:
```bash
conda env create -f environment.yml
conda activate twas-eqtl
```

### R
- `matchit` (>= 4.5.0) — for Mahalanobis matching
- S-PrediXcan (v1.0, MetaXcan v0.8.1) — see [https://github.com/hakyimlab/MetaXcan](https://github.com/hakyimlab/MetaXcan)

## Usage

```bash
# Step 1: Run enrichment analysis and statistical tests
python scripts/python/03_enrichment_analysis.py

# Step 2: Early-generation downstream outputs
python scripts/python/01_visualize_mahalanobis_love_plot.py
python scripts/python/02_density_scatter_consistency.py

# Step 3: The manuscript's current figure values
cd figure_scripts_officialZ_20260917 && python paths_config.py
```

## Reproducibility

The authoritative processed data tables are in `data/processed_officialZ/` (Z-scores, arm denominators, the
four-arm SCZ matrix); the earlier `data/processed/` layer is retained only as pre-correction provenance.
Analysis scripts are version-controlled in this repository under the MIT license.

> ⚠️ `data/processed/eqtlgen_vs_gtex_comparison.csv` is **superseded for direction-consistency analyses**: it
> predates the three-way allele harmonisation, so its `Same_Direction` column reflects pre-harmonisation
> Z-scores. Use `data/processed_officialZ/eqtlgen_official_Z.csv` instead.

The sqrt(N_e)-weighted sensitivity re-merge of the RNH1 DR cross-cohort meta-analysis (M6(b)) and the
descriptive three-study merge re-including GCST90043640 (M6(d)) are reproducible via
`scripts/python/m6_ne_weighted_sensitivity.py` (output: `data/processed/m6_ne_weighted_sensitivity_results.txt`);
the per-gene eQTLGen Z, BH q and FDR calls underlying the group rows of Tables 1b and 2 are provided in
`data/processed_officialZ/eqtlgen_official_Z.csv` (`BH_q` and `FDR_significant`; 96 genes × 3 phenotypes = 288 rows).

## Citation

Wu Y, Chen M, Wu Q, Zhao J, Jin G. (2026). eQTL weight-source choice reshapes TWAS gene candidacy: a two-axis
dual-source audit with disease-agnostic calibration and a genome-wide benchmark (code & data) (v1.0.0).
Zenodo.

<!-- A DOI is minted by Zenodo when the v1.0.0 release is published. The DOI badge above and the DOI in
     this citation are added in the commit that immediately follows that first publication. They are
     deliberately absent here so that no unresolved DOI string ships inside the released archive. -->

## License

MIT License. See [LICENSE](LICENSE) for details.
