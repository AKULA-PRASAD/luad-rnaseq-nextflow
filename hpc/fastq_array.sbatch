#!/bin/bash
# Convert downloaded SRA runs to gzipped paired FASTQ, one array task per run (24 runs).
# Submit from the project folder:  sbatch ~/luad-rnaseq-nextflow/hpc/fastq_array.sbatch
#SBATCH --job-name=luad_fastq
#SBATCH --partition=short
#SBATCH --array=1-24
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/fastq_%a.out
set -euo pipefail
PIPE=${PIPE:-$HOME/luad-rnaseq-nextflow}
SHEET=${SHEET:-$PIPE/assets/samplesheet_GSE37764.csv}
module load sratoolkit/12Dec2024
RUN=$(awk -F, -v n=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==1{for(i=1;i<=NF;i++) if($i=="run") c=i} NR==n{print $c}' "$SHEET")
mkdir -p fastq tmp
if [ -s "fastq/${RUN}_1.fastq.gz" ] && [ -s "fastq/${RUN}_2.fastq.gz" ]; then echo "$RUN already done"; exit 0; fi
fasterq-dump --split-files --threads 4 --temp tmp -O fastq "sra/$RUN"
gzip -f "fastq/${RUN}_1.fastq" "fastq/${RUN}_2.fastq"
ls -la fastq/${RUN}_*
