#!/bin/bash
# =============================================================================
# run_spredixcan.sh — S-PrediXcan 运行脚本（参考/可选）
# =============================================================================
# 此脚本记录了 S-PrediXcan 的完整命令行调用，供完整管道复现使用。
# 需要挂载外部大文件（FinnGen GWAS、GTEx 模型库、1000G LD 参考面板）。
#
# 使用方式:
#   docker run --rm \
#     -v /path/to/finngen_data:/input/finngen \
#     -v /path/to/gtex_models:/input/gtex \
#     -v /path/to/1000g:/input/1000g \
#     -v $(pwd)/output:/app/output \
#     twas-eqtl-repro \
#     bash /app/run_spredixcan.sh
# =============================================================================

set -euo pipefail

INPUT_DIR="/app/input"
OUTPUT_DIR="/app/output"
SPREDIXCAN="/opt/MetaXcan/software/SPrediXcan.py"
COVARIANCE="${INPUT_DIR}/1000g/g1000_eur_covariance.txt.gz"

echo "====================================================================="
echo " S-PrediXcan Execution Script (Reference)"
echo "====================================================================="
echo " Input directory: ${INPUT_DIR}"
echo " Output directory: ${OUTPUT_DIR}"
echo "====================================================================="

# =============================================================================
# PART A: GTEx v8 TWAS (Baseline)
# =============================================================================
echo ""
echo "=== PART A: GTEx v8 TWAS ==="

for TISSUE in "Nerve_Tibial" "Whole_Blood"; do
    MODEL="${INPUT_DIR}/gtex/GTEx_v8_MASHR_${TISSUE}.db"
    for PHENO in "DR" "DN" "DPN"; do
        GWAS="${INPUT_DIR}/finngen/finngen_R13_DM_${PHENO}.gz"
        OUTPUT="${OUTPUT_DIR}/gtex_${TISSUE}_${PHENO}.csv"

        echo "  Running: ${TISSUE} × ${PHENO}..."
        if [ -f "${MODEL}" ] && [ -f "${GWAS}" ]; then
            python3 "${SPREDIXCAN}" \
                --model_db_path "${MODEL}" \
                --covariance "${COVARIANCE}" \
                --gwas_file "${GWAS}" \
                --snp_column rsid \
                --effect_allele_column alt \
                --non_effect_allele_column ref \
                --beta_column beta \
                --se_column se \
                --pvalue_column pval \
                --output_file "${OUTPUT}" \
                --additional_output \
                2>&1 | tail -5
            echo "    Output: ${OUTPUT}"
        else
            echo "    SKIP: missing model or GWAS file"
        fi
    done
done

# =============================================================================
# Multi-tissue integration (Stouffer)
# =============================================================================
echo ""
echo "=== Multi-tissue Stouffer integration ==="
echo "  (Manual step: requires per-tissue output files from above)"
echo "  Formula: Z_combined = sum(w_i * Z_i) / sqrt(sum(w_i^2))"
echo "    where w_i = sqrt(N_i), N_i = tissue sample size"
echo "  Tissues: Nerve_Tibial (N=532) + Whole_Blood (N=670)"

# =============================================================================
# PART B: eQTLGen TWAS (Sensitivity)
# =============================================================================
echo ""
echo "=== PART B: eQTLGen TWAS ==="

MODEL="${INPUT_DIR}/eqtlgen/eQTLGen_Whole_Blood.db"
for PHENO in "DR" "DN" "DPN"; do
    GWAS="${INPUT_DIR}/finngen/finngen_R13_DM_${PHENO}.gz"
    OUTPUT="${OUTPUT_DIR}/eqtlgen_Whole_Blood_${PHENO}.csv"

    echo "  Running: eQTLGen × ${PHENO}..."
    if [ -f "${MODEL}" ] && [ -f "${GWAS}" ]; then
        python3 "${SPREDIXCAN}" \
            --model_db_path "${MODEL}" \
            --covariance "${COVARIANCE}" \
            --gwas_file "${GWAS}" \
            --snp_column rsid \
            --effect_allele_column alt \
            --non_effect_allele_column ref \
            --beta_column beta \
            --se_column se \
            --pvalue_column pval \
            --output_file "${OUTPUT}" \
            --additional_output \
            2>&1 | tail -5
        echo "    Output: ${OUTPUT}"
    else
        echo "    SKIP: missing model or GWAS file"
    fi
done

echo ""
echo "====================================================================="
echo " S-PrediXcan runs complete."
echo "====================================================================="
echo " Next steps: Run downstream analysis with:"
echo "   bash /app/run_all.sh"
echo "====================================================================="
