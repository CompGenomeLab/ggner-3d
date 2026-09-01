#!/bin/bash
#SBATCH --job-name=rnaseq_nf
#SBATCH --output=/cta/users/vkaya/rnaseq/logs/nf_rnaseq_%j.out
#SBATCH --error=/cta/users/vkaya/rnaseq/logs/nf_rnaseq_%j.err
#SBATCH --account=investor
#SBATCH --qos=long_investor
#SBATCH --partition=long_investor
#SBATCH --cpus-per-task=4

cd /cta/users/vkaya/rnaseq
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
# Run nextflow
NXF_OPTS='-Xms1g -Xmx4g -Duser.language=en -Duser.region=US' nextflow run /cta/users/vkaya/rnaseq/rnaseq/main.nf \
  --input ./spreadsheet.csv \
  --outdir ./results \
  --fasta /cta/users/vkaya/chipseq_work/refs/GRCh38.primary_assembly.genome.fa \
  --gtf /cta/users/vkaya/chipseq_work/refs/gencode.v48.primary_assembly.annotation.gtf \
  -profile slurm,conda
