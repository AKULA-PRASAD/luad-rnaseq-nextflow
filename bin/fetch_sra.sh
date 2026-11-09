#!/usr/bin/env bash
# Download the SRA runs listed in a samplesheet (column `run`).
# Run on a node WITH internet access (on many clusters, the login node).
# Usage: bin/fetch_sra.sh assets/samplesheet_GSE37764.csv <outdir>
set -euo pipefail
sheet=${1:?samplesheet}; out=${2:?output directory}
mkdir -p "$out/sra"
awk -F, 'NR==1{for(i=1;i<=NF;i++) if($i=="run") c=i; next} {print $c}' "$sheet" | while read -r run; do
  prefetch "$run" -O "$out/sra" --max-size 60G
done
echo "SRA files in $out/sra; convert with hpc/fastq_array.sbatch (or fasterq-dump --split-files)"
