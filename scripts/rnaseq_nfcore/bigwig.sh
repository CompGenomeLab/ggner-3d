#!/usr/bin/env bash
#SBATCH --job-name=rnaseq_merge_bw
#SBATCH --output=/cta/users/vkaya/rnaseq/logs/rnaseq_merge_bw_%A_%a.out
#SBATCH --error=/cta/users/vkaya/rnaseq/logs/rnaseq_merge_bw_%A_%a.err
#SBATCH --account=adelab
#SBATCH --qos=adelab
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=32
#SBATCH --array=0-3

set -euo pipefail

ROOT_DIR="/cta/users/vkaya/rnaseq"
STAR_DIR="${ROOT_DIR}/results/star_salmon"
OUT_DIR="${ROOT_DIR}/results/merged_tracks"

MAMBA_ENV="${MAMBA_ENV:-deeptools}"
ENV_PREFIX="${ENV_PREFIX:-/cta/users/vkaya/micromamba/envs/${MAMBA_ENV}}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
BIN_SIZE="${BIN_SIZE:-100}"
NORMALIZATION="${NORMALIZATION:-CPM}"

DEFAULT_CONDITIONS=(RNA_0 RNA_12 RNA_30 RNA_60)

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#DEFAULT_CONDITIONS[@]} )); then
        echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is out of range." >&2
        exit 1
    fi
    CONDITIONS=( "${DEFAULT_CONDITIONS[${SLURM_ARRAY_TASK_ID}]}" )
elif [[ $# -gt 0 ]]; then
    CONDITIONS=("$@")
else
    echo "Submit with sbatch ${BASH_SOURCE[0]} or pass one or more condition names." >&2
    exit 1
fi

SAMTOOLS_BIN="${SAMTOOLS_BIN:-${ENV_PREFIX}/bin/samtools}"
BAMCOVERAGE_BIN="${BAMCOVERAGE_BIN:-${ENV_PREFIX}/bin/bamCoverage}"

if [[ ! -x "${SAMTOOLS_BIN}" ]]; then
    echo "samtools not found at ${SAMTOOLS_BIN}" >&2
    exit 1
fi

if [[ ! -x "${BAMCOVERAGE_BIN}" ]]; then
    echo "bamCoverage not found at ${BAMCOVERAGE_BIN}" >&2
    exit 1
fi

"${SAMTOOLS_BIN}" --version >/dev/null
"${BAMCOVERAGE_BIN}" --version >/dev/null

mkdir -p "${OUT_DIR}/bam" "${OUT_DIR}/bigwig"

for condition in "${CONDITIONS[@]}"; do

    merged_bam="${OUT_DIR}/bam/${condition}.merged.markdup.sorted.bam"
    merged_bw="${OUT_DIR}/bigwig/${condition}.merged.${NORMALIZATION,,}.v2.bigWig"
    # Exclude unmapped, secondary, QC-fail, and supplementary alignments.
    "${BAMCOVERAGE_BIN}" \
        --bam "${merged_bam}" \
        --outFileName "${merged_bw}" \
        --binSize "${BIN_SIZE}" \
        --normalizeUsing "${NORMALIZATION}" \
        --numberOfProcessors "${THREADS}"
done

printf 'Normalized bigWigs: %s/bigwig\n' "${OUT_DIR}"