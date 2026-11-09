#!/usr/bin/env nextflow
/*
 * LUAD RNA-seq: raw FASTQ -> QC -> STAR alignment -> gene counts -> DESeq2
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
    if (!(params.strandedness in ['auto', 'unstranded', 'forward', 'reverse']))
        error "--strandedness must be auto, unstranded, forward or reverse"

    def sheet = file(params.samplesheet, checkIfExists: true)
    def base  = sheet.parent

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

    MULTIQC(
        FASTQC.out.zip.map { _meta, z -> z }
            .mix(STAR_ALIGN.out.log.map { _meta, l -> l })
            .collect()
    )
}
