#!/bin/bash
# Quick check on the cluster: bundled test data, containers from the local cache.
#SBATCH --job-name=luad_nf_test
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/nextflow_test.out
set -euo pipefail
PIPE=${PIPE:-$HOME/lung-cancer-rnaseq-nextflow}
module load nextflow/25.04.6
export NXF_SINGULARITY_CACHEDIR=$PWD/containers NXF_OFFLINE=true NXF_ANSI_LOG=false
mkdir -p test_run && cd test_run
nextflow run "$PIPE/main.nf" -profile test,singularity
python3 "$PIPE/tests/check_outputs.py" results_test
