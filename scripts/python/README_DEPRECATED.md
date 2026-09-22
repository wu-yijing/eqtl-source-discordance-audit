# ⚠️ EARLY-GENERATION pipeline — do not use to reproduce the paper

The scripts in this directory belong to the **early (iScience-era) generation** of this project.
They are retained for provenance only.

## Why they are not the current pipeline

- They read `data/processed/` — the **pre-correction** intermediate layer (produced before the
  S-PrediXcan σᵢ expression-variance factor and the PLINK 2-bit decoding defect were fixed).
  The authoritative layer is `data/processed_officialZ/`.
- They write their figures to `figs/`, which is a **runtime** directory (created by `run_all.sh`
  and then copied to `output/`), not part of the archived figure set.
- Their numbers therefore disagree with the published figures — for example the eQTLGen-arm
  Z-scores of the candidate genes (RNH1 13.318 vs 2.3064; TUBB 48.516 vs 11.8932;
  CKAP4 4.7549 vs 0.9874; and the GTEx-side RNH1 value +13.82 vs +2.67).

## What to use instead

| Purpose | Use |
|---|---|
| The manuscript's figure values | `figure_scripts_officialZ_20260917/` (start with `python paths_config.py`); output goes to a `.gitignore`d `figures/` |
| Z-scores, denominators, cross-cohort values | `data/processed_officialZ/` (with `_PROVENANCE.json`) |
| Genome-wide SCZ arm | Superseded in-house snapshot — **not distributed with this repository** (the current SCZ statistics come from the official MetaXcan v0.8.1 binary; see `data/processed_officialZ/scz_z_4arm_official.csv`) |

`run_all.sh` still invokes `scripts/python/04_generate_all_figures.py`; that container path
reproduces the earlier pipeline and the environment, not the current figure set.

*(Added 2026-09-20 as part of the layout cleanup.)*
