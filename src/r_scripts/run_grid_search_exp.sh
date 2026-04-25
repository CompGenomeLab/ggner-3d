#!/bin/bash

# --- Configuration ---
# Set global constants
BASE_COVARIATES=""
CRITICAL_COLS=""
R_SCRIPT="./xpc_ctcf_grid_search_3.R"
PROJECT_ROOT="/home/carlos/Clone/ggner-3d"
FOCAL_KEPT_FRACT=(1.00 0.90 0.80 0.70 0.60 0.50 0.40 0.30 0.20 0.10)
# FOCAL_KEPT_FRACT=(1.00)
OUT_FOLDER="1_0_1"
GRID_WORKERS="${GRID_WORKERS:-0}"

# --- Function Definition ---
run_grid_search() {
    local ins_ws="$1"
    local extra_covs="$2"

    # Combine covariates
    local covariates="${BASE_COVARIATES} ${extra_covs}"
    
    # Create a safe filename suffix (replace spaces with underscores)
    local suffix="${extra_covs// /_}"

    # Define paths
    local input_file="${PROJECT_ROOT}/src/r_scripts/WTnoUV_ws_${ins_ws}_with_pol2.csv"
    local output_dir="${PROJECT_ROOT}/data/${OUT_FOLDER}/ws_${ins_ws}"
    
    local summary_file="${output_dir}/grid_search_summary_${suffix}.csv"
    local smd_file="${output_dir}/grid_search_smds_${suffix}.csv"
    local output_data="${output_dir}/grid_search_data_${suffix}.csv"

    # Ensure directories exist
    mkdir -p "$output_dir"

    echo "---------------------------------------------------"
    echo "Running analysis for WS: ${ins_ws} | Extra: ${extra_covs}"
    echo "Output: ${output_dir}"
    echo "Grid workers: ${GRID_WORKERS}"
    
    # Run the R script
    Rscript "$R_SCRIPT" \
      --input_file $input_file \
      --summary_file $summary_file \
      --smd_file $smd_file \
      --covariates $covariates \
      --critical_exact_match_cols $CRITICAL_COLS \
      --drop_failed_strata \
      --output_data $output_data \
      --workers $GRID_WORKERS \
      --focal_keep_frac_seq "${FOCAL_KEPT_FRACT[@]}"
}

# --- Execution ---

run_grid_search "100000" "pol2"
run_grid_search "100000" "atac_signal"
run_grid_search "100000" "h3k4me3_signal"
run_grid_search "100000" "h3k27ac_signal"
run_grid_search "100000" "h3k4me1_signal"
