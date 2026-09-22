# S1 cluster-robustness analysis — eQTL-source discordance in TWAS

> ⚠️ **Superseded generation (2026-09-20).** Every number on this page was computed on the
> **pre-correction** Z layer (`data/processed/`) and describes the v2.5.0-generation run — e.g.
> 61/96 = 63.5% and ρ = 0.321, where the current manuscript has 66/96 = 68.8% and ρ = 0.39.
> The cluster-aware procedures the manuscript now reports (naive t = 4.10, df 94, P = 8.8 × 10⁻⁵;
> sandwich SE = 0.125; jackknife SE = 0.137; permutation P < 0.001; gene-cluster bootstrap ρ CI
> 0.12–0.62 and direction consistency 58.3–79.2%) are in **Additional file 1: Table S16** and in the
> repository README section *Cluster-Robustness Re-analysis of Direction Consistency (S1)*.
> This directory is retained for provenance only; the code below is not the script set behind those
> numbers and should not be used to reproduce them.

Code and results for the gene-level **cluster-robust re-analysis of direction consistency**
reported in the manuscript *"eQTL-source discordance in TWAS: a dual-source sensitivity
analysis with true-negative calibration"*.

Purpose: the primary between-source comparison comprises 96 gene–phenotype pairs that
derive from only 32 genes. Because a gene's three phenotypes (DR, DN, DPN) are not
independent observations, the originally reported two-sided exact binomial test
(P = 0.010) could in principle be anti-conservative. This directory contains the
code that (i) reconstructs the analysis arms exactly, (ii) quantifies the clustering,
and (iii) recomputes the test with three cluster-aware procedures.

---

## 1. Data provenance

All inputs are the analysed output tables already archived in this Zenodo record:

| File | Role |
|---|---|
| `data/processed/gtex_Whole_Blood_{DR,DN,DPN}.csv` | GTEx v8 MASHR whole-blood S-PrediXcan Z |
| `data/processed/gtex_Nerve_Tibial_{DR,DN,DPN}.csv` | GTEx v8 MASHR tibial-nerve Z |
| `data/processed/eqtlgen_spredixcan_harmonized_results.csv` | eQTLGen Z after three-way allele harmonisation (**use this file**) |
| `data/processed/covariate_matrix.csv` | gene-group assignments (candidate / non-candidate / T2DM control) |
| `_DEPRECATED_scz_self_implemented/results/scz_twas_results_limit0.csv` | genome-wide SCZ arm: 10,357 genes × 4 Z columns — ⚠️ **not distributed with this repository**; supply your own copy via `S1_SCZ_SNAPSHOT` |
| `_DEPRECATED_scz_self_implemented/results/scz_axis_difference.json` | reference values for the SCZ decomposition — ⚠️ **not distributed** |

> ⚠️ **Two inputs are no longer shipped.** The in-house SCZ snapshot (`n = 2,511`) was retired when
> the SCZ statistics were recomputed with the official MetaXcan v0.8.1 binary, so the SCZ part of
> this directory cannot be re-run from the repository alone — set `S1_SCZ_SNAPSHOT` (see
> `s1_scz_cluster.py`) to your own copy if you need it. The path names above are kept as written
> at the time for audit traceability.

> ⚠️ **Do not use `data/processed/eqtlgen_vs_gtex_comparison.csv`** for this analysis.
> That file (dated 2026-09-02) predates the three-way allele harmonisation and its
> `Same_Direction` column reflects the pre-harmonisation eQTLGen Z. Recomputing the
> 34-gene anchor set from it gives 61/102 (59.8%) instead of the reported 64/102 (62.7%).
> All numbers in this directory were produced from the harmonised file.
> The stale file is retained in the record for provenance only.

---

## 2. Reconstruction of the analysis arms

**Primary arm (32 genes / 96 pairs)** — non-candidate (44) and T2DM-control (30) genes
for which both sources yield a Z in all three phenotypes; genes with a missing
eQTLGen Z are dropped. Of the 36 eligible genes, four have `NaN` eQTLGen Z across all
three phenotypes (GCKR, HSPD1, KLF14, PCBP2 — 2 non-candidates + 2 T2DM controls),
leaving **21 non-candidates + 11 T2DM controls = 32 genes / 96 pairs**. This matches
the arm described in the manuscript.

**Anchor set (34 genes / 102 pairs)** — the primary arm minus NCL and PDIA6, plus the
four non-testbed proteins (SERPINH1, RPL13, VAT1, RPL17).

**Full testbed (52 genes / 156 pairs)** — all groups (candidate + non-candidate + T2DM
control) with valid Z from both sources, dropping genes with missing eQTLGen Z.

### Reproduction check (scripts vs. values in the manuscript)

| Arm | Reproduced here | Manuscript | Match |
|---|---|---|---|
| Primary | 61/96 = 63.54%; ρ = 0.321 (P = 0.0014); binomial P = 0.0103 | 61 (63.5%); ρ = 0.32 (P = 0.001); P = 0.010 | ✅ |
| Primary, per phenotype | DR 21/32, DN 21/32, DPN 19/32 | 21/32, 21/32, 19/32 | ✅ |
| Anchor set | 64/102 = 62.7%; ρ = 0.266 (P = 0.0068); P = 0.0129 | 64 (62.7%); ρ = 0.27 (P = 0.007); P = 0.013 | ✅ |
| Full testbed | 109/156 = 69.9% | 109 (69.9%) | ✅ |
| SCZ source axis | 1,714/2,511 = 68.26%; ρ = 0.5222 | 68.3%; ρ = 0.52 | ✅ |
| SCZ tissue axis | 1,788/2,511 = 71.21%; ρ = 0.5088 | 68.3% / 71.2% | ✅ |

---

## 3. Statistical methods

1. **Design effect / intra-cluster correlation (ICC).** For the binary consistency
   indicator, a one-way ANOVA decomposition over gene clusters:
   `ICC = (MSB − MSW) / [MSB + (m̄ − 1)·MSW]`, `m̄` = 3; `DEFF = 1 + (m̄ − 1)·ICC`.
2. **Cluster-robust (sandwich) variance** of the ratio estimator
   `p = Σk_g / Σn_g`:
   `Var(p) = [K/(K−1)]·Σ(k_g − p·n_g)² / (Σn_g)²`, `K` = number of genes;
   `z = (p − 0.5)/√Var`, two-sided normal P.
3. **Gene-level cluster bootstrap.** Genes resampled with replacement (B = 10,000);
   all of a resampled gene's pairs enter the resample. Percentile 95% CI and
   `P = 2 × min[P(p* ≤ 0.5), P(p* ≥ 0.5)]`.
4. **Gene-label permutation test.** Within each phenotype, eQTLGen Z-scores are
   permuted across genes (B = 10,000), preserving marginal distributions and per-gene
   pair counts. This tests H0: no cross-gene correspondence between the two sources.
5. **Spearman ρ under gene resampling.** Same gene-level bootstrap, reporting the
   percentile 95% CI of ρ.

SCZ arm note: the complete-case set contributes **exactly one pair per gene**
(2,511 pairs from 2,511 genes), so the design effect is 1 by construction and the
gene-resampling bootstrap coincides with the exact binomial test.

---

## 4. Results

### 4.1 HOTAIR testbed arms

| Arm | genes / pairs | consistency | exact binomial P | ICC | DEFF | cluster-robust z P | gene-cluster bootstrap P (95% CI) | permutation P | ρ (bootstrap 95% CI) |
|---|---|---|---|---|---|---|---|---|---|
| Primary | 32 / 96 | 61 (63.5%) | 0.0103 | 0.000 | 1.000 | **0.0133** | **0.0168** (53.1–74.0%) | **0.0078** | 0.321 (0.039–0.561) |
| Anchor set | 34 / 102 | 64 (62.7%) | 0.0129 | 0.000 | 1.000 | 0.0226 | 0.0322 (52.0–73.5%) | 0.0064 | 0.266 (−0.010–0.516) |
| Full testbed | 52 / 156 | 109 (69.9%) | <1e−4 | 0.021 | 1.043 | <1e−4 | <1e−4 (60.9–78.2%) | <1e−4 | 0.447 (0.237–0.622) |

**Conclusion.** The intra-gene correlation is negligible (ICC = 0.000, DEFF = 1.000),
so the exact binomial test reported in the manuscript is **not materially
anti-conservative**: the cluster-robust and permutation P-values span 0.008–0.017
against the reported 0.010, and the 95% CI widens by only 0.9 percentage points
(53.1–74.0% vs 53.1–73.1%). All variants remain significant at α = 0.05.

### 4.2 Genome-wide SCZ arm (2,511 pairs)

| Axis | pairs | consistency | exact binomial P | DEFF | bootstrap 95% CI | ρ |
|---|---|---|---|---|---|---|
| resource/sample-size (eQTLGen WB vs GTEx multi-tissue) | 2,511 | 1,714 (68.26%) | 2.4e−76 | 1.0 | 66.43–70.05% | 0.5222 |
| tissue-context (GTEx Whole_Blood vs Nerve_Tibial) | 2,511 | 1,788 (71.21%) | 3.2e−103 | 1.0 | 69.41–72.96% | 0.5088 |

---

## 5. Files

| File | Description |
|---|---|
| `s1_recon.py` | reconstructs the primary arm from raw per-source tables; prints arm composition |
| `s1_diag.py` | diagnoses the NaN genes and the per-gene cluster structure |
| `s1_fix.py` | confirms the four `NaN` eQTLGen genes and rebuilds the 96-pair primary arm |
| `s1_cluster.py` | **main analysis**: ICC/DEFF, cluster-robust variance, gene-cluster bootstrap, gene-label permutation for the three HOTAIR arms |
| `s1_scz_cluster.py` | genome-wide SCZ arm re-analysis (shows DEFF = 1 by construction) |
| `s1_results.json` | machine-readable results, HOTAIR arms |
| `s1_scz_results.json` | machine-readable results, SCZ arm |
| `s1_primary_arm_96pairs.csv` | the 96 primary-arm pairs with Z_GTEx, Z_eQTLGen and Same (direction) — allows independent verification |
| `s1_anchor_102pairs.csv` | the 102 anchor-set pairs |

## 6. Environment and reproducibility

- Python 3.13.12; `numpy`, `pandas`, `scipy` (`scipy.stats.spearmanr`, `binomtest`, `norm`).
  `statsmodels` is not used.
- Random seeds: `numpy.random.default_rng(20260910)` for the HOTAIR arms
  (10,000 resamples); `numpy.random.default_rng(20260726)` for the SCZ arm, matching
  the seed used in `_DEPRECATED_scz_self_implemented/scz_axis_difference_analysis.py`.
- Run order and working directory: place this directory alongside the repository root
  so that the relative `data/processed/…` paths resolve, then

  ```bash
  python s1_recon.py       # verify arm composition
  python s1_diag.py        # inspect NaN genes / cluster sizes
  python s1_cluster.py     # main table (Section 4.1)
  python s1_scz_cluster.py # SCZ table (Section 4.2)
  ```

  Alternatively edit the `D` / `SRC` constants at the top of each script to absolute paths.
- Wall time: `s1_cluster.py` ≈ 3 min (10,000 × 2 bootstrap/permutation loops);
  `s1_scz_cluster.py` ≈ 5 min for 2 × 10,000 iterations over 2,511 observations.
