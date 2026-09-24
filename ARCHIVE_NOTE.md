# ARCHIVE_NOTE — how this archive maps onto the manuscript

**Manuscript:** "eQTL weight-source choice reshapes TWAS gene candidacy: a two-axis dual-source audit with disease-agnostic calibration and a genome-wide benchmark" — submitted to *BMC Genomics*, 24 September 2026.

**Archive version:** v1.0.1. Between v1.0.0 and v1.0.1 **no data file and no computed value changed**; the release corrects documentation only (see *What changed in v1.0.1*). Every number reported in the manuscript is therefore reproducible from **either** version.

**Scope.** This repository is a code-and-processed-data release. It deliberately contains **no figures, figure captions, tables or table legends** — those accompany the manuscript and its Additional file 1. Raw RNA pull-down / LC–MS/MS spectra are deposited separately in ProteomeXchange via iProX.

---

## 1. Which files reproduce which manuscript items

All paths are relative to the repository root.

| Manuscript item | Authoritative file(s) | Notes |
|---|---|---|
| Testbed gene list with group assignments | `data/processed_officialZ/gene_groups_TableS1_official.csv` (104 rows) | = Additional file 1: Table S1 |
| GTEx v8 MASHR multi-tissue per-gene Z / P / FDR q | `data/processed_officialZ/gtex_official_Z.csv` (222 rows = 74 genes × 3 phenotypes) | = Additional file 1: Table S2 |
| Primary 96-pair comparison (both source Z per pair) | `data/processed_officialZ/primary_arm_96pairs_official.csv` (96 rows) | = Additional file 1: Table S12; headline ρ = 0.39 and 68.8% direction consistency reproduce here |
| eQTLGen whole-blood per-gene Z / P / BH q / model SNPs | `data/processed_officialZ/eqtlgen_official_Z.csv` (288 rows, 96 genes) | The 69-gene universe of Additional file 1: Tables S17/S18 is a subset |
| Cross-cohort values (RNH1, FinnGen vs UK Biobank) | `data/processed_officialZ/crosscohort_TableS4_official.csv` | = Additional file 1: Table S4a |
| Genome-wide SCZ four-arm decomposition | `data/processed_officialZ/scz_z_4arm_official.csv` (15,875 rows) | 8,315-gene complete-case set; = Additional file 1: Table S23 |
| Framework-only / allele-alignment contrasts | `data/processed_officialZ/` + `figure_scripts_officialZ_20260917/` | = Additional file 1: Table S15 |
| Four anchor-set protein statistics from the earlier in-house implementation | `data/processed/gtex_acat_o_results.csv` and `gtex_stouffer_integrated.csv` (**237 rows = 79 genes × 3**) | This is the **"archived working table"** named in Additional file 1: Table S10. Superseded layer — see §2 |
| Independence of the matched-vs-pool signal-density check | `data/processed/mahalanobis_matched_pairs.csv` + `data/processed_officialZ/*` | Matched pairs are covariate-only; see §3 |

## 2. Authoritative vs superseded layers

| Path | Status |
|---|---|
| `data/processed_officialZ/` | ✅ **Authoritative** — official MetaXcan v0.8.1 (S-PrediXcan reference implementation) after three-way allele harmonisation |
| `data/processed/` | ⚠️ **Superseded** — pre-correction intermediate layer (missing S-PrediXcan σᵢ expression-variance factor; PLINK 2-bit decoding defect). Record: `data/processed/_PROVENANCE.json` |
| `scripts/python/` | ⚠️ **Early generation** — reads the superseded layer and writes figures to a runtime `figs/`. Record: `scripts/python/README_DEPRECATED.md` |
| `figure_scripts_officialZ_20260917/` | ✅ The only pipeline that reproduces the manuscript's current figure values |
| `data/hk_reselect_20260830/` | ✅ Housekeeping-control layer (Layer 1 and its reselection) |

The two defects in the earlier in-house implementation are disclosed in the manuscript's Methods; the earlier routine is retained only as an equivalence cross-check, where the residual difference from the official binary was at most |ΔZ| = 3 × 10⁻⁸.

## 3. How Additional file 1: Table S25 is reproduced

Table S25 (Mahalanobis-matched enrichment contrasts, added at the 24 Sep 2026 revision) requires three inputs, all of them in this archive:

1. `data/processed/mahalanobis_matched_pairs.csv` — the 30 candidate–control matched pairs (`subclass`, `treated`);
2. `data/processed_officialZ/gtex_official_Z.csv` — GTEx endpoint (`FDR_q_ACAT_O` for the BH endpoint, `P_ACAT_O` for the nominal endpoint);
3. `data/processed_officialZ/eqtlgen_official_Z.csv` — eQTLGen endpoint (`BH_q` and `P`).

Denominators are the testable pairs: GTEx 84 candidate vs 60 matched-control pairs (28 vs 20 genes); eQTLGen 81 vs 57 (27 vs 19 genes). Two-sided Fisher exact tests give 1.00 and 0.077 (BH endpoint, GTEx and eQTLGen) and 0.78 and 0.76 (nominal endpoint).

**Why the matched-pairs file is usable despite being in the superseded layer.** It contains covariates only (gene length, GC content, eQTL SNP counts, matched subclass). The two S-PrediXcan defects act on the *Z-score* computation and cannot propagate to covariates; the 30 pairs and every numeric covariate value are identical to those in Additional file 1: Table S3.

## 4. Container / pipeline entry point — read this before running anything

`run_all.sh` is retained for environment provenance. **Its default path executes the early-generation scripts** in `scripts/python/` (see `scripts/python/README_DEPRECATED.md`), which read `data/processed/` and therefore print pre-correction values — including an enrichment rate for the GTEx candidate arm that does not match the manuscript.

To reproduce manuscript values, use `figure_scripts_officialZ_20260917/` (start with `python paths_config.py`) and read `data/processed_officialZ/`. In particular, `scripts/python/03_enrichment_analysis.py` emits a report whose conclusions predate the implementation correction (for example, it states that the matched candidate-versus-control Fisher test could not be performed — it can, and the result is Additional file 1: Table S25).

## 5. What changed in v1.0.1 (documentation only)

| File | Change |
|---|---|
| `data/README.md` | Rewritten: superseded manuscript title removed, GigaDB framing removed, **licence corrected from CC0-1.0 to MIT** to match the manuscript and the Zenodo record, row counts corrected to the archived values, dangling pointer to a non-distributed `analysis_reports/` removed, authoritative-vs-superseded tables added |
| `ARCHIVE_NOTE.md` | Added (this file) |
| `README.md` | Pointer to this file added under "Repository scope" |
| `scripts/python/README_DEPRECATED.md` | Named `03_enrichment_analysis.py` explicitly and pointed at this file |
| `.zenodo.json` | `version` 1.0.0 → 1.0.1 |

**No file under `data/` was modified**, and no computed value reported in the manuscript changed between v1.0.0 and v1.0.1.

## 6. Licence

Code **and** processed data: **MIT** (see `LICENSE`), matching the manuscript's Data availability statement.
