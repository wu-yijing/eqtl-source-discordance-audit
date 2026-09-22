# HRT Atlas v1.0 source gene universe (archived for reproducibility)

`Human_Mouse_Common.csv` is the byte-faithful source file of the housekeeping-gene universe
used for the disease-agnostic control (Methods 1.8) and for the architecture-unselected
genome-wide/HRT controls (Additional file 1: Tables S10-S11).

* Original URL : https://housekeeping.unicamp.br/Human_Mouse_Common.csv
* Reference    : Hounkpe BW, Chenou F, de Lima F, De Paula EV. HRT Atlas v1.0 database:
  redefining human and mouse housekeeping genes and candidate reference transcripts by
  mining massive RNA-seq datasets. Nucleic Acids Res. 2021;49(D1):D947-D955. PMID 32663312
* Mirror       : MSigDB gene set HOUNKPE_HOUSEKEEPING_GENES (systematic name M42508)
* Format       : semicolon-delimited `Mouse;Human`; 1,130 data rows -> 1,129 unique human symbols

## Verification fingerprint (reproduced 2026-09-11)

Running `scripts/01_select_hk_genes.py` rules R1-R5 on this file reproduces the recorded
filter chain exactly:

| step | description | n |
|---|---|---|
| R1 | HRT Atlas human-mouse common HK set | 1129 |
| R2 | minus 104-gene testbed panel overlap | 1112 |
| R3 | minus families sharing a >=3-char leading prefix with a panel gene (+ MRPS*/MRPL*/MT-*) | 1013 |
| R4 | minus curated T2DM / diabetic-complication / insulin & GLP-1 / metabolic genes | 1007 |
| R5 | require a MASHR model with >=1 SNP in BOTH Nerve_Tibial and Whole_Blood | 767 |

All 30 genes of `hk_genes_v2.txt` are contained in the R5 pool.

## Second control (Additional file 1: Table S11)

Relaxing R5 to ">=1 SNP in Whole_Blood only" (no second-tissue requirement) gives the
architecture-unselected HRT pool of **818 genes** used for the HRT-restricted random control.
