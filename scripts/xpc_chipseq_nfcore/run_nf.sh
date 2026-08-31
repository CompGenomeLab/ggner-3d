#!/bin/bash
#SBATCH --job-name=xpc_chipseq
#SBATCH --output=/cta/users/vkaya/chipseq_work/xpc_results/logs/nf_chipseq_%j.out
#SBATCH --error=/cta/users/vkaya/chipseq_work/xpc_results/logs/nf_chipseq_%j.err
#SBATCH --account=investor
#SBATCH --qos=long_investor
#SBATCH --partition=long_investor
#SBATCH --cpus-per-task=4

set -euo pipefail

cd /cta/users/vkaya/chipseq_work
eval "$("${HOME}/miniconda3/bin/conda" shell.bash hook)"

NXF_OPTS='-Xms1g -Xmx4g -Duser.language=en -Duser.region=US' \
  "${HOME}/.local/bin/nextflow" -trace nextflow.executor \
  run chipseq-2.1.0/main.nf \
  --input ./xpc_spreadsheet.csv \
  --outdir ./xpc_results \
  --genome GRCh38 \
  --save_reference \
  -profile slurm,conda \
  --narrow_peak \
  --macs_gsize 2913022398 \
  --save_macs_pileup \
  -resume
