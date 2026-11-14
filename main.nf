#!/usr/bin/env nextflow
/*
 * LUAD RNA-seq: raw FASTQ -> QC -> STAR alignment -> gene counts -> DESeq2 -> GSEA
 *
 * Lung adenocarcinoma tumour vs matched normal tissue from never-smoker women
 * (GEO GSE37764 / SRA SRP012656: 6 patients, paired-end RNA-seq, 2 runs per sample).
 */
nextflow.enable.dsl = 2

include { FASTQC       } from './modules/fastqc'
include { STAR_INDEX   } from './modules/star_index'
include { STAR_ALIGN   } from './modules/star_align'
include { MERGE_COUNTS } from './modules/merge_counts'
include { DESEQ2       } from './modules/deseq2'
include { GSEA         } from './modules/gsea'
include { MULTIQC      } from './modules/multiqc'

// Resolve a samplesheet path: absolute paths and URLs as given, relative paths
// against the folder that holds the samplesheet.
def resolvePath(String p, base) {
    if (!p) return null
    return (p.startsWith('/') || p.contains('://')) ? file(p, checkIfExists: true)
                                                   : file(base.resolve(p).toString(), checkIfExists: true)
}

workflow {
    if (!params.samplesheet) error "Set --samplesheet (columns: sample,condition,fastq_1[,fastq_2,run,<batch>])"
    if (!params.star_index && !(params.genome && params.gtf)) error "Set --genome and --gtf, or --star_index"
    if (!params.gtf) error "Set --gtf (needed for gene names and counting)"
    if (!params.skip_gsea && !params.gsea_gmt) error "Set --gsea_gmt, or --skip_gsea"
    if (!(params.strandedness in ['auto', 'unstranded', 'forward', 'reverse']))
        error "--strandedness must be auto, unstranded, forward or reverse"

    def sheet = file(params.samplesheet, checkIfExists: true)
    def base  = sheet.parent

    // Whole-sheet checks that a per-row check cannot see
    def rows = sheet.splitCsv(header: true, strip: true)
    if (!rows) error "Samplesheet ${sheet} has no rows"
    def fastqs = rows.collectMany { r -> [r.fastq_1, r.fastq_2].findAll { v -> v } }
    def dup_fq = fastqs.countBy { v -> v }.findAll { _k, n -> n > 1 }.keySet()
    if (dup_fq) error "Samplesheet lists the same FASTQ file more than once: ${dup_fq.join(', ')}"
    def runs = rows.collect { r -> r.run ?: r.fastq_1 }
    if (runs.unique(false).size() != runs.size()) error "Samplesheet has duplicate values in the run column"

    // One row per sequencing run. Several rows with the same `sample` are merged at alignment.
    ch_rows = Channel.fromPath(sheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.condition || !row.fastq_1)
                error "Samplesheet row missing sample, condition or fastq_1: ${row}"
            if (params.batch_column && !row[params.batch_column])
                error "Samplesheet row has no value in batch column '${params.batch_column}': ${row}"
            def meta = [id: row.sample, condition: row.condition,
                        batch: params.batch_column ? row[params.batch_column] : '']
            def run  = row.run ?: row.sample
            def r1   = resolvePath(row.fastq_1, base)
            def r2   = row.fastq_2 ? resolvePath(row.fastq_2, base) : null
            tuple(row.sample, meta, run, r1, r2)
        }

    // FastQC per run
    ch_fastqc = ch_rows.map { _s, _meta, run, r1, r2 -> tuple([id: run], r2 ? [r1, r2] : [r1]) }

    // STAR per sample: all runs of a sample together, read 1 and read 2 kept apart
    ch_samples = ch_rows
        .map { s, meta, _run, r1, r2 -> tuple(s, meta, r1, r2 ?: 'NONE') }
        .groupTuple()
        .map { s, metas, r1s, r2s ->
            if (metas.unique().size() != 1) error "Sample ${s}: rows disagree on condition or batch"
            def paired = r2s.every { it != 'NONE' }
            if (!paired && r2s.any { it != 'NONE' }) error "Sample ${s}: mixes single-end and paired-end runs"
            tuple(metas[0] + [single_end: !paired], r1s, paired ? r2s : [])
        }

    gtf = file(params.gtf, checkIfExists: true)
    index = params.star_index
        ? Channel.value(file(params.star_index, checkIfExists: true))
        : STAR_INDEX(file(params.genome, checkIfExists: true), gtf).index

    FASTQC(ch_fastqc)
    STAR_ALIGN(ch_samples, index)

    MERGE_COUNTS(
        STAR_ALIGN.out.counts.map { _meta, tab -> tab }.collect(),
        gtf,
        params.strandedness
    )

    DESEQ2(MERGE_COUNTS.out.matrix, sheet)

    if (!params.skip_gsea) {
        GSEA(DESEQ2.out.results, file(params.gsea_gmt, checkIfExists: true))
    }

    MULTIQC(
        FASTQC.out.zip.map { _meta, z -> z }
            .mix(STAR_ALIGN.out.log.map { _meta, l -> l })
            .collect()
    )

    // One software_versions.yml for the run: each process reports its tool versions
    ch_versions = FASTQC.out.versions
        .mix(STAR_ALIGN.out.versions, MERGE_COUNTS.out.versions, DESEQ2.out.versions, MULTIQC.out.versions)
    if (!params.star_index) ch_versions = ch_versions.mix(STAR_INDEX.out.versions)
    if (!params.skip_gsea)  ch_versions = ch_versions.mix(GSEA.out.versions)
    ch_versions
        .map { f -> f.text }
        .unique()
        .mix(Channel.of("\"Workflow\":\n    luad-rnaseq-nextflow: \"${workflow.manifest.version}\"\n    nextflow: \"${workflow.nextflow.version}\"\n".toString()))
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)
}
