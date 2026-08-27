#!/bin/bash
#SBATCH --job-name=chipseq_main
#SBATCH --output=/cta/users/vkaya/chipseq_work/ctcf/logs/nf_chipseq_%j.out
#SBATCH --error=/cta/users/vkaya/chipseq_work/ctcf/logs/nf_chipseq_%j.err
#SBATCH --account=investor
#SBATCH --qos=long_investor
#SBATCH --partition=long_investor
#SBATCH --cpus-per-task=4

cd /cta/users/vkaya/chipseq_work/ctcf
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
fasta="/cta/users/vkaya/chipseq_work/refs/GRCh38.primary_assembly.genome.fa"
gtf="/cta/users/vkaya/chipseq_work/refs/gencode.v48.primary_assembly.annotation.gtf"
aligner="chromap"
index="/cta/users/vkaya/chipseq_work/refs/chromap_index/GRCh38.primary_assembly.genome.index"
# Run nextflow
NXF_OPTS='-Xms1g -Xmx4g -Duser.language=en -Duser.region=US' \
  $HOME/.local/bin/nextflow run /cta/users/vkaya/chipseq_work/chipseq-2.1.0/main.nf \
  --input ./spreadsheet.csv \
  --outdir ./results \
  --fasta $fasta \
  --gtf $gtf \
  --aligner $aligner \
  --chromap_index $index \
  --save_reference \
  -profile slurm,conda \
  --narrow_peak \
  --macs_gsize 2913022398 \
  --save_macs_pileup \
  -resume
  # -trace nextflow.executor  # before run
