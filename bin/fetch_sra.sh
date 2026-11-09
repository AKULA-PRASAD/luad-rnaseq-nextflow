#!/usr/bin/env bash
# Download SRA runs and convert to gzipped FASTQ.
# Run on a node WITH internet access (on many HPC clusters, the login node),
# because compute nodes often cannot reach NCBI.
#
# Usage: bin/fetch_sra.sh assets/samplesheet_SRP012656.csv <outdir>
# Needs SRA Toolkit (prefetch, fasterq-dump) on PATH, e.g. `module load sratoolkit`
# or `conda install -c bioconda sra-tools=3.2.0`.
set -euo pipefail
sheet=${1:?samplesheet}; out=${2:?output directory}
mkdir -p "$out/fastq" "$out/sra"
tail -n +2 "$sheet" | cut -d, -f1 | while read -r run; do
  [ -s "$out/fastq/$run.fastq.gz" ] && { echo "skip $run (exists)"; continue; }
  prefetch "$run" -O "$out/sra"
  fasterq-dump --split-files --threads 4 -O "$out/fastq" "$out/sra/$run"
  # single-end runs give $run.fastq; paired-end runs give _1/_2 files
  for f in "$out/fastq/$run"*.fastq; do gzip -f "$f"; done
done
echo "FASTQ files in $out/fastq"
