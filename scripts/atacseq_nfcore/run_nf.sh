#!/bin/bash
#SBATCH --job-name=atac_nf
#SBATCH --output=/cta/users/vkaya/atac/logs/nf_atac_%j.out
#SBATCH --error=/cta/users/vkaya/atac/logs/nf_atac_%j.err
#SBATCH --account=investor
#SBATCH --qos=long_investor
#SBATCH --partition=long_investor
#SBATCH --cpus-per-task=4

cd /cta/users/vkaya/atac
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
fasta="/cta/users/vkaya/chipseq_work/refs/GRCh38.primary_assembly.genome.fa"
gtf="/cta/users/vkaya/chipseq_work/refs/gencode.v48.primary_assembly.annotation.gtf"
aligner="chromap"
index="/cta/users/vkaya/chipseq_work/refs/chromap_index/GRCh38.primary_assembly.genome.index"
# Run nextflow
NXF_OPTS='-Xms1g -Xmx4g -Duser.language=en -Duser.region=US' \
nextflow run /cta/users/vkaya/atac/atacseq/main.nf \
  --input ./spreadsheet_xpc.csv \
  --outdir ./results \
  --fasta "$fasta" \
  --gtf "$gtf" \
  --aligner "$aligner" \
  --chromap_index "$index" \
  -profile slurm,conda \
  --narrow_peak \
  --macs_gsize 2913022398 \
  --save_macs_pileup \
  --blacklist /cta/users/vkaya/atac/hg38-blacklist.v2.bed \
  -resume
  # -trace nextflow.executor  # before run
