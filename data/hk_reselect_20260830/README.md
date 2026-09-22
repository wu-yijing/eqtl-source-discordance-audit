# Housekeeping true-negative control set — v2 provenance

This directory archives the full provenance of the **v2 housekeeping (HK) true-negative
control gene set** used in the manuscript (Table S6, Methods §*, and Figure 9).

## Why v2 exists

The v1 set (`data/hk_genes.txt` as committed before 2026-08-30) contained **7 of 30 genes
that overlapped the 104-gene testbed panel** — `ACTB`, `GAPDH`, `RPL13A`, `RPLP0`, `RPS18`
(candidates) and `EEF2`, `ENO1` (non-candidates) — which directly violated its own
selection criterion #3 ("Not present in the study's existing 104-gene panel"). v2 removes
all of them and re-samples from a clean eligible pool.

## v2 selection rules (pre-specified, outcome-blind)

- **R1 universe** — HRT Atlas v1.0 human–mouse common HK set (MSigDB
  `HOUNKPE_HOUSEKEEPING_GENES`, n = 1129). Ref: Hounkpe BW et al., *Nucleic Acids Res.*
  2021;49(D1):D947–D955. PMID 32663312.
- **R2 exclude** — any gene in the 104-gene testbed panel.
- **R3 exclude** — any gene sharing a ≥3-char leading alphabetic prefix with a 104-panel
  gene (family-level exclusion: `HSP*`, `RPS*`, `RPL*`, `EEF*`, `ENO*`, `ANXA*`, `ATP*`,
  `HNRNP*`, `PDIA*`, `ILF*`, `XRCC*`, `SLC*`, `YWHA*`, …), plus `MRPS*`/`MRPL*`/`MT-*`.
- **R4 exclude** — curated T2DM / diabetic-complication / insulin- & GLP-1-signalling genes.
- **R5 require** — GTEx v8 MASHR model with ≥1 SNP in **both** `Nerve_Tibial` and
  `Whole_Blood`.
- **R6 sampling** — simple random sample of 30, `random.seed(20260830)`.

### Filter flow
| Step | Remaining |
|------|-----------|
| R1 HRT Atlas universe | 1129 |
| R2 minus 104-panel overlap | 1112 |
| R3 minus same-family (≥3-char prefix) | 1013 |
| R4 minus T2DM/complication/metabolic | 1007 |
| R5 require dual-tissue MASHR model | 767 |
| R6 random sample (seed 20260830) | 30 |

## Testability note (28 of 30)

All 30 genes are listed in `data/hk_genes.txt` and Table S6. Of these, **28 are GTEx-testable**
(i.e. have ≥1 imputed TWAS Z in both `Nerve_Tibial` and `Whole_Blood` under the 1000G EUR
LD reference). The two exceptions — **`SPRYD3`** and **`TUT1`** — lack usable LD-panel
variants in the GTEx v8 MASHR models for one or both tissues, so they cannot yield a TWAS
Z and are reported as "not testable" in Table S6. The pairing analyses therefore use n = 28
HK pairs (see `hk_v2_summary.json` → `enrichment_*` → `N_tested` = 28 for the Housekeeping v2 rows).

## Key outputs

- `data/hk_genes_v2.txt` — canonical 30-gene list (identical to `data/hk_genes.txt`).
- `TableS6_hk_control_v2.csv` — Table S6 source data (30 genes, testability + SNP counts).
- `Table2_v2.csv` — Table 2 HK-v2 enrichment cells.
- `hk_selection_meta.json` — seed, pool sizes, full gene list, per-gene model SNP counts.
- `hk_v2_summary.json` — FDR-enrichment tables (ACAT-O / Nerve_Tibial / Whole_Blood) and
  Figure 9 statistics.
- `arms_all_groups.csv` — per-gene TWAS arms used by Figure 9 panel c.
- `panel104.txt` / `validate_panel104.csv` — 104-panel reference + overlap validation.
- `hk_twas_v2_raw.csv` — raw S-PrediXcan HK-v2 results.

### Figure 9 (v2)
- panel (a) GTEx multi-tissue (ACAT-O) vs Nerve_Tibial: Spearman ρ = 0.95 (95% CI 0.92–0.97), n = 75
- panel (b) GTEx Whole_Blood vs Nerve_Tibial: Spearman ρ = 0.68 (95% CI 0.52–0.79), n = 66
- panel (c) tissue-context axis across gene sets: PGC SCZ3 +0.51 (n = 2511) / HOTAIR testbed
  +0.45 (n = 150) / Housekeeping v2 +0.68 (n = 66)

## Reproducing the scripts

The `scripts/` folder contains the exact code run on 2026-08-30. **These scripts have
hard-coded local Windows paths** (`ROOT`, `REPO`, `HRT`, `MODEL_DIR`) and depend on data
not shipped with this repo (the HRT Atlas CSV and GTEx v8 MASHR model directory). To re-run,
edit those constants at the top of each script to point at your environment. Steps:

1. `01_select_hk_genes.py` — applies R1–R6 and writes `hk_genes_v2.txt`.
2. `02_spredixcan.py` — runs S-PrediXcan on the v2 set.
3. `03_enrichment_and_figs.py` — enrichment + supplementary figures.
4. `04_regenerate_table2.py` — Table 2 HK-v2 cells.
5. `05_ld_panel_audit.py` — 1000G EUR LD-panel audit (identifies SPRYD3/TUT1 as not testable).
6. `06_figure9_v2.py` — Figure 9 (three panels).
7. `verify_all.py` — end-to-end consistency checks.
