#!/usr/bin/env Rscript
# =============================================================================
# run_mahalanobis_matching.R
# Mahalanobis covariate matching for TWAS candidate gene enrichment validation
# =============================================================================
# This script reproduces the matching of 30 HOTAIR candidate genes against
# a background gene pool using MatchIt.
#
# Input:  data/processed/covariate_matrix.csv  (pool of genes with covariates)
# Output: data/processed/mahalanobis_matched_pairs.csv  (30 matched pairs)
#         figs/Fig5_Love_Plot.pdf  (covariate balance diagnostic)
# =============================================================================

# ---- 0. Setup ----
SCRIPT_DIR <- dirname(normalizePath(sys.frame(1)$ofile))
PROJECT_DIR <- file.path(SCRIPT_DIR, "..", "..")
DATA_DIR <- file.path(PROJECT_DIR, "data", "processed")
FIG_DIR <- file.path(PROJECT_DIR, "figs")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load data ----
cat("[Mahalanobis] Loading covariate matrix...\n")
covar <- read.csv(file.path(DATA_DIR, "covariate_matrix.csv"),
                  stringsAsFactors = FALSE)

# =============================================================================
# NOTE (2026-09-11): provenance of missing eQTL SNP counts
# -----------------------------------------------------------------------------
# The archived outputs (mahalanobis_matched_pairs.csv, Additional file 1
# Table S4) were generated with group-median imputation of eQTL_SNPs_Mean for
# genes with no significant cis-eQTL model in the corresponding weight-source
# run: 3 of 30 candidates (HSPA8, RPS18, RPS25; median 1.5), 17 of 44 pool
# genes (10 of which entered the matched control set; median 1.5), and 11 of
# 30 T2DM control genes (median 1.0). The na.omit-style filtering in Section 2
# below therefore does NOT reproduce those outputs as-is; it is retained to
# document the complete-case sensitivity view. To reproduce the published
# matched pairs, impute eQTL_SNPs_Mean with the median of observed values
# within each Group (Candidate / Non-Candidate / T2DM Control) BEFORE the
# filtering steps. See manuscript Methods 1.9 for the full disclosure.
# =============================================================================

# ---- 2. Prepare matching data ----
# Select candidate genes (treated)
candidate <- covar[covar$Group == "Candidate", ]
candidate <- candidate[!is.na(candidate$eQTL_SNPs_Mean), ]

# Create matching pool: all non-candidate + T2DM control genes
# Exclude candidate genes from the pool
pool <- covar[covar$Group != "Candidate", ]
pool <- pool[!is.na(pool$eQTL_SNPs_Mean), ]

# Prepare matching dataframe
# Add treated flag: 1 = candidate, 0 = background
matching_data <- rbind(
    data.frame(
        Gene = candidate$Gene,
        Group = "Treated",
        log10_Length = log10(candidate$Length_bp),
        GC_pct = candidate$GC_pct,
        eQTL_SNPs_Mean = candidate$eQTL_SNPs_Mean,
        treated = 1,
        stringsAsFactors = FALSE
    ),
    data.frame(
        Gene = pool$Gene,
        Group = "Control",
        log10_Length = log10(pool$Length_bp),
        GC_pct = pool$GC_pct,
        eQTL_SNPs_Mean = pool$eQTL_SNPs_Mean,
        treated = 0,
        stringsAsFactors = FALSE
    )
)

cat(sprintf("  Treated (candidates): %d\n", sum(matching_data$treated == 1)))
cat(sprintf("  Control pool:         %d\n", sum(matching_data$treated == 0)))

# Remove rows with missing values
matching_data <- na.omit(matching_data)
cat(sprintf("  After NA removal:     %d treated, %d controls\n",
    sum(matching_data$treated == 1),
    sum(matching_data$treated == 0)))

# ---- 3. Run Mahalanobis matching ----
cat("[Mahalanobis] Running MatchIt (Mahalanobis, nearest neighbor, ratio=1)...\n")

library(MatchIt)
library(cobalt)

m.out <- matchit(treated ~ log10_Length + GC_pct + eQTL_SNPs_Mean,
                 data = matching_data,
                 method = "nearest",
                 distance = "mahalanobis",
                 ratio = 1,
                 replace = FALSE)

m.data <- match.data(m.out)

cat(sprintf("  Matched pairs: %d / %d (%.0f%%)\n",
    nrow(m.data[m.data$treated == 1, ]),
    sum(matching_data$treated == 1),
    nrow(m.data[m.data$treated == 1, ]) / sum(matching_data$treated == 1) * 100))

# ---- 4. Extract matched pairs ----
pairs <- m.data[, c("Gene", "Group", "log10_Length", "GC_pct",
                     "eQTL_SNPs_Mean", "subclass", "treated")]
colnames(pairs)[colnames(pairs) == "log10_Length"] <- "log10_Length"
colnames(pairs)[colnames(pairs) == "eQTL_SNPs_Mean"] <- "n_eQTL_SNPs"
# Add Length_bp back
pairs$Length_bp <- round(10^pairs$log10_Length)

# Write output
output_file <- file.path(DATA_DIR, "mahalanobis_matched_pairs.csv")
write.csv(pairs, output_file, row.names = FALSE)
cat(sprintf("  Output written: %s\n", output_file))

# ---- 5. Balance diagnostics (Love Plot) ----
cat("[Mahalanobis] Generating Love Plot (covariate balance diagnostic)...\n")

pdf(file.path(FIG_DIR, "Fig5_Love_Plot.pdf"), width = 8, height = 5)
love.plot(m.out,
          binary = "std",
          stats = c("mean.diffs"),
          threshold = 0.1,
          abs = TRUE,
          var.order = "unadjusted",
          line = TRUE,
          colors = c("#E74C3C", "#3498DB"),
          shapes = c(18, 16),
          sample.names = c("Before Matching", "After Matching"),
          title = "Covariate Balance: Mahalanobis Matching",
          labels = c(
              log10_Length = "log10(Gene Length)",
              GC_pct = "GC Content (%)",
              eQTL_SNPs_Mean = "eQTL SNPs (mean)"
          ))
dev.off()
cat("  Love Plot saved: figs/Fig5_Love_Plot.pdf\n")

# ---- 6. Print SMD summary ----
cat("\n[Mahalanobis] Standardized Mean Differences:\n")
summary_matched <- summary(m.out, standardize = TRUE)
print(summary_matched$sum.all[, c("Diff.Adj", "Diff.Unadj")])

cat("\n[Mahalanobis] Matching complete.\n")
