#!/bin/bash
#SBATCH --account=mdbf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --output=/cta/users/vkaya/chipseq_work/xpc_results/logs/bw_norm_%A_%a.out
#SBATCH --error=/cta/users/vkaya/chipseq_work/xpc_results/logs/bw_norm_%A_%a.err
#SBATCH --job-name=xpc_bw_norm
#SBATCH --array=0-12               # adjust to match number of BAMs

set -euo pipefail

# --- USER EDITS --------------------------------------------------------------
BASE="/cta/users/vkaya/chipseq_work/xpc_results"
BAMS="${BASE}/bwa/merged_library"
OUT="${BASE}/bw_deeptools"
mkdir -p "${OUT}/rpgc" "${OUT}/ratio" "${BASE}/logs"

# Genome / RPGC
EGS=2913022398                      # hg38 effective genome size
BLACKLIST="${BASE}/hg38-blacklist.v2.bed"

# Coverage params
BINSIZE=25
SMOOTH=75
FRAGLEN=180
THREADS=${SLURM_CPUS_PER_TASK}

# Tools env
source ~/micromamba/etc/profile.d/mamba.sh
micromamba activate deeptools     # contains deepTools + samtools + macs3

# Build the ordered BAM list
mapfile -t BAMS_ARR < <(ls -1 ${BAMS}/XPC_*sorted.bam | sort)
N=${#BAMS_ARR[@]}
IDX=${SLURM_ARRAY_TASK_ID}
echo "Index ${IDX} of ${N} in ${BAMS}"
if [[ ${IDX} -ge ${N} ]]; then
  echo "Index ${IDX} >= ${N}; nothing to do"; exit 0
fi

BAM="${BAMS_ARR[$IDX]}"
BN=$(basename "${BAM%.sorted.bam}")

# Index BAM if needed
[[ -f "${BAM}.bai" ]] || samtools index -@ ${THREADS} "${BAM}"

# --- Detect PE vs SE.
PAIRED_READS=$(samtools view -c -f 1 "${BAM}" || echo 0)
if [[ "${PAIRED_READS}" -gt 0 ]]; then
  IS_PE=1
else
  IS_PE=0
fi

# --- 1) Per-sample RPGC bigWig
COV_CMD=( bamCoverage
  -b "${BAM}"
  -o "${OUT}/rpgc/${BN}.rpgc.bw"
  --normalizeUsing RPGC
  --effectiveGenomeSize ${EGS}
  --binSize ${BINSIZE} --smoothLength ${SMOOTH}
  --blackListFileName "${BLACKLIST}"
  --ignoreDuplicates
  --numberOfProcessors ${THREADS}
)

# Only extend for SE
if [[ "${IS_PE}" -eq 0 ]]; then
  COV_CMD+=( --extendReads "${FRAGLEN}" )
fi

echo "[${BN}] Running: ${COV_CMD[*]}"
"${COV_CMD[@]}"

# --- 2) Per-replicate WT vs matching KO (SES log2 ratio)
# Infer timepoint by filename tokens: 1h / 3h / noUV
if [[ "${BN}" == XPC_WT_* ]]; then
  if   [[ "${BN}" == *"_1h_"* ]];   then KO="${BAMS}/XPC_KO_1h_REP1.mLb.clN.sorted.bam"   ; LAB="1h"
  elif [[ "${BN}" == *"_3h_"* ]];   then KO="${BAMS}/XPC_KO_3h_REP1.mLb.clN.sorted.bam"   ; LAB="3h"
  elif [[ "${BN}" == *"_noUV_"* ]]; then KO="${BAMS}/XPC_KO_noUV_REP1.mLb.clN.sorted.bam"; LAB="noUV"
  else
    echo "[${BN}] Could not infer timepoint; skipping ratio."
    exit 0
  fi

  [[ -f "${KO}.bai" ]] || samtools index -@ ${THREADS} "${KO}"

  bamCompare \
    -b1 "${BAM}" -b2 "${KO}" \
    -o "${OUT}/ratio/${BN}_vs_KO_${LAB}.log2.SES.bw" \
    --operation log2 --pseudocount 1 \
    --scaleFactorsMethod SES \
    --binSize ${BINSIZE} --smoothLength ${SMOOTH} \
    --blackListFileName "${BLACKLIST}" \
    --ignoreDuplicates \
    --numberOfProcessors ${THREADS}
fi
