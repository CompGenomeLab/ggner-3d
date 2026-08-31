#!/bin/bash
#SBATCH --account=adelab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --qos=adelab
#SBATCH --partition=genomics
#SBATCH -w cn20
#SBATCH --job-name=atac_bw_rpgc
#SBATCH --array=0-1
#SBATCH --output=/cta/users/vkaya/atac/logs/bw_rpgc_%A_%a.out
#SBATCH --error=/cta/users/vkaya/atac/logs/bw_rpgc_%A_%a.err

set -euo pipefail

source ~/micromamba/etc/profile.d/mamba.sh
micromamba activate deeptools

# ===============================
# PATHS / CONSTANTS (hg38)
# ===============================
BASEDIR="/cta/users/vkaya/atac/results/chromap/merged_replicate"
BLACKLIST="/cta/users/vkaya/atac/hg38-blacklist.v2.bed"
EFFECTIVE_GENOME_SIZE=2913022398   # hg38 deepTools effective genome size

# ===============================
# DEFAULT PARAMETERS (override via args)
# Usage:
#   sbatch atac_bw_rpgc.slurm [binSize] [smoothLength] [minMAPQ] [ATACshift 0/1]
# Example:
#   sbatch atac_bw_rpgc.slurm 10 30 30 1
# ===============================
BIN_SIZE=10
SMOOTH_LEN=30
MIN_MAPQ=30
DO_ATACSHIFT=0

if [[ $# -ge 1 ]]; then BIN_SIZE=$1; fi
if [[ $# -ge 2 ]]; then SMOOTH_LEN=$2; fi
if [[ $# -ge 3 ]]; then MIN_MAPQ=$3; fi
if [[ $# -ge 4 ]]; then DO_ATACSHIFT=$4; fi

echo "Parameters:"
echo "  binSize=${BIN_SIZE}"
echo "  smoothLength=${SMOOTH_LEN}"
echo "  minMAPQ=${MIN_MAPQ}"
echo "  ATACshift=${DO_ATACSHIFT}"

# ===============================
# INPUTS (array)
# ===============================
# BAMS=(
#   "${BASEDIR}/WT_noUV.mRp.clN.sorted.bam"
#   "${BASEDIR}/WT_UV3h.mRp.clN.sorted.bam"
# )

BAMS=(
  "${BASEDIR}/XPC_noUV.mRp.clN.sorted.bam"
  "${BASEDIR}/XPC_UV3h.mRp.clN.sorted.bam"
)

BAM="${BAMS[${SLURM_ARRAY_TASK_ID}]}"
SAMPLE="$(basename "${BAM%.bam}")"

OUTDIR="${BASEDIR}/bigwig/rpgc"
TMPDIR="${BASEDIR}/bigwig/tmp_shift"
mkdir -p "${OUTDIR}" "${TMPDIR}"

# ===============================
# Optional ATAC Tn5 shift (must re-sort!)
# ===============================
BAM_FOR_COV="${BAM}"

if [[ "${DO_ATACSHIFT}" -eq 1 ]]; then
  SHIFTED_BAM="${TMPDIR}/${SAMPLE}.ATACshift.unsorted.bam"
  SHIFTED_SORTED_BAM="${TMPDIR}/${SAMPLE}.ATACshift.sorted.bam"

  if [[ ! -s "${SHIFTED_SORTED_BAM}" ]]; then
    echo "[${SAMPLE}] alignmentSieve --ATACshift..."
    alignmentSieve \
      -b "${BAM}" \
      --ATACshift \
      -p "${SLURM_CPUS_PER_TASK}" \
      -o "${SHIFTED_BAM}"

    echo "[${SAMPLE}] samtools sort shifted BAM..."
    samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "${SHIFTED_SORTED_BAM}" "${SHIFTED_BAM}"

    echo "[${SAMPLE}] samtools index shifted+sorted BAM..."
    samtools index -@ "${SLURM_CPUS_PER_TASK}" "${SHIFTED_SORTED_BAM}"

    rm -f "${SHIFTED_BAM}"
  fi

  BAM_FOR_COV="${SHIFTED_SORTED_BAM}"
fi

# ===============================
# bamCoverage (RPGC normalized bigWig)
# ===============================
OUTBW="${OUTDIR}/${SAMPLE}.RPGC.bs${BIN_SIZE}.sl${SMOOTH_LEN}.mapq${MIN_MAPQ}.ATACshift${DO_ATACSHIFT}.bw"

echo "[${SAMPLE}] bamCoverage -> ${OUTBW}"
bamCoverage \
  -b "${BAM_FOR_COV}" \
  -o "${OUTBW}" \
  -p "${SLURM_CPUS_PER_TASK}" \
  --binSize "${BIN_SIZE}" \
  --smoothLength "${SMOOTH_LEN}" \
  --minMappingQuality "${MIN_MAPQ}" \
  --extendReads \
  --normalizeUsing RPGC \
  --effectiveGenomeSize "${EFFECTIVE_GENOME_SIZE}" \
  --blackListFileName "${BLACKLIST}" \
  --ignoreForNormalization chrX chrY chrM

echo "[${SAMPLE}] Done."
echo "  Input BAM: ${BAM_FOR_COV}"
echo "  Output BW: ${OUTBW}"
