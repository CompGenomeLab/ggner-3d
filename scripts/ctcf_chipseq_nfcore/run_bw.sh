#!/bin/bash
#SBATCH --account=investor
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --qos=mid_investor
#SBATCH --partition=mid_investor
#SBATCH --output=/cta/users/vkaya/chipseq_work/ctcf/logs/bw_norm_%A_%a.out
#SBATCH --error=/cta/users/vkaya/chipseq_work/ctcf/logs/bw_norm_%A_%a.err
#SBATCH --job-name=ctcf_bw_norm
#SBATCH --array=0-8                    # 9 BAMs -> 0..8 (use %N to throttle if needed, e.g. 0-8%3)

set -euo pipefail

# --- USER EDITS --------------------------------------------------------------
BASE="/cta/users/vkaya/chipseq_work/ctcf/results"
BAMS="${BASE}/chromap/merged_library"
OUT="${BASE}/bw_deeptools"
mkdir -p "${OUT}/rpgc" "${OUT}/ratio" "${BASE}/logs"

# Genome / normalization (hg38)
EGS=2913022398
BLACKLIST="/cta/users/vkaya/chipseq_work/ctcf/hg38-blacklist.v2.bed" # Boyle Lab blacklist (hg38)

# Coverage parameters (CTCF = narrow TF)
BINSIZE=25
SMOOTH=75
THREADS=${SLURM_CPUS_PER_TASK}

# Env
source ~/micromamba/etc/profile.d/mamba.sh && micromamba activate deeptools
# ------------------------------------------------------------------------------

# Ordered BAM list (exact files you showed, then sorted for determinism)
mapfile -t BAMS_ARR < <(ls -1 \
  ${BAMS}/INPUT_REP1.mLb.clN.sorted.bam \
  ${BAMS}/WT_3h_REP1.mLb.clN.sorted.bam \
  ${BAMS}/WT_3h_REP2.mLb.clN.sorted.bam \
  ${BAMS}/WT_NoUV_REP1.mLb.clN.sorted.bam \
  ${BAMS}/WT_NoUV_REP2.mLb.clN.sorted.bam \
  ${BAMS}/XPCKO_3h_REP1.mLb.clN.sorted.bam \
  ${BAMS}/XPCKO_3h_REP2.mLb.clN.sorted.bam \
  ${BAMS}/XPCKO_NoUV_REP1.mLb.clN.sorted.bam \
  ${BAMS}/XPCKO_NoUV_REP2.mLb.clN.sorted.bam | sort)

N=${#BAMS_ARR[@]}
IDX=${SLURM_ARRAY_TASK_ID}
echo "Index ${IDX} of ${N} in ${BAMS}"
if [[ ${IDX} -ge ${N} ]]; then
  echo "Index ${IDX} >= ${N}; nothing to do"; exit 0
fi

BAM="${BAMS_ARR[$IDX]}"
BN=$(basename "${BAM%.sorted.bam}")

# Shared Input
INPUT="${BAMS}/INPUT_REP1.mLb.clN.sorted.bam"

# Index BAMs if needed
[[ -f "${BAM}.bai" ]]   || samtools index -@ ${THREADS} "${BAM}"
[[ -f "${INPUT}.bai" ]] || samtools index -@ ${THREADS} "${INPUT}"

# 1) Per-sample RPGC bigWig (PE => no --extendReads)
bamCoverage \
  -b "${BAM}" \
  -o "${OUT}/rpgc/${BN}.rpgc.bw" \
  --normalizeUsing RPGC \
  --effectiveGenomeSize ${EGS} \
  --binSize ${BINSIZE} --smoothLength ${SMOOTH} \
  --blackListFileName "${BLACKLIST}" \
  --ignoreDuplicates \
  --numberOfProcessors ${THREADS}

# 2) Per-replicate ChIP vs shared Input ratio (SES log2)
#    Skip if current BAM is the Input itself
if [[ "${BN}" != INPUT_* ]]; then
  bamCompare \
    -b1 "${BAM}" -b2 "${INPUT}" \
    -o "${OUT}/ratio/${BN}_vs_INPUT.log2.SES.bw" \
    --operation log2 --pseudocount 1 \
    --scaleFactorsMethod SES \
    --binSize ${BINSIZE} --smoothLength ${SMOOTH} \
    --blackListFileName "${BLACKLIST}" \
    --ignoreDuplicates \
    --numberOfProcessors ${THREADS}
fi
