# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-11-14

First public release.

### Added
- Nextflow DSL2 pipeline with these steps: FastQC, STAR index, STAR alignment with gene counts, count merging, DESeq2 and MultiQC.
- Pre-ranked gene set enrichment (fgsea) on the DESeq2 Wald statistic, with the MSigDB Hallmark gene sets bundled (`--gsea_gmt`, `--skip_gsea`).
- Samplesheet with one row per sequencing run. Runs that share a sample are aligned together, and paired-end and single-end data are both supported.
- Optional batch term in the DESeq2 design (`--batch_column`), for example `~ patient + condition` for matched tumour and normal tissue.
- Automatic strandedness detection from STAR's per-strand gene counts.
- Every process records its tool versions, collected in `pipeline_info/software_versions.yml`.
- Validation of the whole samplesheet: the run stops if a FASTQ file or run is listed twice, if a sample's rows disagree on condition or batch, or if single-end and paired-end runs are mixed.
- Profiles for Docker, Singularity, Apptainer, conda and SLURM, and two for Northeastern's Explorer cluster (`explorer`, `explorer_node`).
- Synthetic test datasets with planted effects (single-end, and paired-end with several runs per sample), checked against the known truth.
- Unit tests for count merging and strandedness detection.
- GitHub Actions CI on Nextflow 24.10.5 and 25.04.8. It runs lint, structure and syntax checks and the unit tests, checks that every profile resolves, tests the samplesheet validation, runs both test profiles, checks the test data is reproducible and checks that `-resume` uses the cache.
- Full analysis of GEO GSE37764 (6 patients, tumour and matched normal), with results and the MultiQC report in `docs/`.

### Fixed during development
- The STAR biocontainer crashed with `Illegal instruction` on some cluster CPUs. The Explorer profiles therefore use the cluster's STAR module.
- Parallel FastQC tasks crashed with `SIGBUS` when they shared `/tmp`. The FastQC step now turns off the Java monitoring file.
- With one job per step, intermittent SLURM submission failures stopped runs. The single-allocation profile `explorer_node` avoids them.

### Data provenance
- An earlier, non-Nextflow version of this analysis used SRA runs SRR502127–SRR502138, which are ENCODE ChIP-seq, not this study's RNA-seq. Its conclusions are withdrawn. See the README.
