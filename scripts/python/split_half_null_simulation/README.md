# Split-half null simulation for the 96-pair primary arm

This directory contains the implementation and results of the split-half null
simulation reported in Methods 2.5 and Results 3.3 of the manuscript
("the operative null ... (ii) by split-half simulation on the primary arm itself").

## What it answers

The cross-source direction-consistency rate (63.5%, 61/96) cannot be tested
against a 50% chance level, because both eQTL weight sources are evaluated
against the same GWAS summary statistics. The split-half simulation supplies the
complementary internal benchmark: **how reproducible is a single S-PrediXcan sign
within one weight source?**

For each of the 32 primary-arm genes, the eQTLGen cis-SNP set is randomly
partitioned into two disjoint halves; the TWAS statistic is recomputed on each
half from the same GWAS and the same LD reference. The two estimates therefore
share the GWAS input and the LD reference but use disjoint SNPs, so their sign
concordance measures estimator reproducibility rather than source disagreement.

## Result

| Stratum | Sign concordance | 95% interval across partitions |
|---|---|---|
| Overall (96 pairs) | **96.69%** | 93.75 – 98.96% |
| DR | 93.84% | 87.50 – 100% |
| DN | 99.06% | 96.80 – 100% |
| DPN | 97.17% | 93.75 – 100% |

B = 200 random partitions, seed 20260910 (`numpy.random.default_rng`).

Interpretation: within a single eQTL weight source the S-PrediXcan sign is
reproducible in ~97% of split-half comparisons, whereas the cross-source rate is
63.5%. The ~36% cross-source reversal rate therefore reflects the weight source
rather than estimation noise.

## Inputs

| Input | Source used here |
|---|---|
| eQTL weights | eQTLGen whole-blood cis-eQTL catalogue, FDR < 0.05 (`AssessedAllele`, `Zscore` columns) |
| GWAS Z | FinnGen R13 DR / DN / DPN, `Z = beta / sebeta` |
| LD reference | 1000 Genomes Phase 3 European panel, N = 503 (PLINK bed/bim/fam) |
| Gene list | the 32 primary-arm genes (21 non-candidates + 11 T2DM controls) |

## Estimator

`Z = (w' z) / sqrt(w' R w)`

* `w` = eQTL Z-score, sign-aligned to the LD reference allele
* `z` = GWAS Z-score, sign-aligned to the same allele
* `R` = SNP–SNP LD correlation matrix from the 1000G EUR panel

For a SNP half `S`, the quadratic form is computed as `u' R u` with
`u = w * 1[S]`, which is algebraically identical to `w_S' R_S w_S` and avoids
materialising sub-matrices.

Allele harmonisation is three-way (eQTLGen assessed/other, FinnGen ref/alt,
1000G A1/A2), with strand-flip handling; unresolvable (palindromic or
mismatched) SNPs are dropped per gene–phenotype pair.

## Reproducing

`split_half_null_simulation.py` reads from absolute local paths by default
(weights, GWAS, PLINK reference). Edit the constants at the top of the file
(`BASE`, `BED`, `BIM_MAP`, `CIS`, `_gwas_*.csv`) to point at your local copies
before running. The GWAS tables are expected as
`rsid, chr, pos, ref, alt, beta, sebeta, Z` and the weights as
`gene, snp, chr, pos, assessed, other, zscore, FDR`.

Intermediate artefacts produced here are not committed (they exceed repository
size limits); `splithalf_per_pair.csv` and `splithalf_summary.csv` are the
archived outputs, and `splithalf_result_raw.json` holds the per-partition
boolean vectors for all 96 pairs.

## Honest caveat

The per-pair eQTLGen Z-scores archived in
`scripts/python/s1_cluster_robustness/s1_primary_arm_96pairs.csv` (the values
used for the 63.5% cross-source rate) could **not** be reproduced with this
pipeline. Thirteen combinations of weight transform (Z-score, ±/× per-SNP SD,
eQTL beta and variants) and estimator form (S-PrediXcan `w'z/√(w'Rw)`, the
repository's GLS `w'R⁻¹z/√(w'R⁻¹w)`, and LD-free `w'z/√(w'w)`) reproduced the
sign in every case but differed in magnitude by a gene-specific constant factor.
The split-half concordance reported here is invariant to such per-gene scaling,
so it remains valid as an internal reproducibility benchmark, but it is not
derived from byte-identical reproduction of the archived per-pair Z-scores.
