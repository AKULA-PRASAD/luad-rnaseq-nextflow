#!/bin/bash
# Full GSE37764 run (6 patients, tumour + matched normal, 24 paired-end runs).
# This job runs the Nextflow head process; each pipeline step is submitted by
# Nextflow as its own SLURM job (profile explorer).
#SBATCH --job-name=luad_nf
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/nextflow.out
set -euo pipefail
PIPE=${PIPE:-$HOME/luad-rnaseq-nextflow}
module load nextflow/25.04.6
export NXF_SINGULARITY_CACHEDIR=$PWD/containers NXF_OFFLINE=true NXF_ANSI_LOG=false
cp "$PIPE/assets/samplesheet_GSE37764.csv" .
# Reuse a STAR index in ./star_index if present (built with sjdbOverhang 77), otherwise build it
if [ -d star_index ]; then REF="--star_index star_index"; else REF="--genome ref/GRCh38.primary_assembly.genome.fa.gz"; fi
nextflow run "$PIPE/main.nf" -profile explorer \
    --samplesheet samplesheet_GSE37764.csv \
    --batch_column patient \
    $REF \
    --gtf ref/gencode.v44.primary_assembly.annotation.gtf.gz \
    --sjdb_overhang 77 \
    --outdir results -resume
