#!/bin/bash
# Convert downloaded SRA runs to gzipped FASTQ, one array task per sample.
# Submit from the project folder:  sbatch ~/luad-rnaseq-nextflow/hpc/fastq_array.sbatch
#SBATCH --job-name=luad_fastq
#SBATCH --partition=short
#SBATCH --array=1-12
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/fastq_%a.out
set -euo pipefail
PIPE=${PIPE:-$HOME/luad-rnaseq-nextflow}
module load sratoolkit/12Dec2024
RUN=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$PIPE/assets/samplesheet_SRP012656.csv" | cut -d, -f1)
mkdir -p fastq tmp
if [ -s "fastq/$RUN.fastq.gz" ]; then echo "$RUN already done"; exit 0; fi
fasterq-dump --split-files --threads 4 --temp tmp -O fastq "sra/$RUN"
# single-end runs may be written as RUN.fastq or RUN_1.fastq; normalise to RUN.fastq
[ -f "fastq/${RUN}_1.fastq" ] && [ ! -f "fastq/${RUN}_2.fastq" ] && mv "fastq/${RUN}_1.fastq" "fastq/$RUN.fastq"
gzip -f fastq/${RUN}*.fastq
ls -la fastq/${RUN}*
