#!/bin/bash
# Full SRP012656 run. This job runs the Nextflow head process; each pipeline step
# is submitted by Nextflow as its own SLURM job (profile slurm,singularity).
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
cp "$PIPE/assets/samplesheet_SRP012656.csv" .
nextflow run "$PIPE/main.nf" -profile slurm,singularity \
    --samplesheet samplesheet_SRP012656.csv \
    --genome ref/GRCh38.primary_assembly.genome.fa.gz \
    --gtf ref/gencode.v44.primary_assembly.annotation.gtf.gz \
    --star_simd sse4.1 \
    --outdir results -resume
