#!/usr/bin/env nextflow
/*
 * LUAD RNA-seq: raw FASTQ -> QC -> STAR alignment -> gene counts -> DESeq2
 *
 * Tumour vs adjacent-normal lung adenocarcinoma (SRA study SRP012656, 12 samples),
 * rebuilt as a Nextflow DSL2 pipeline from the parameters of the original HPC run.
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
    if (!params.samplesheet) error "Set --samplesheet (columns: sample,condition,fastq_1[,fastq_2])"
    if (!params.star_index && !(params.genome && params.gtf)) error "Set --genome and --gtf, or --star_index"
    if (!params.gtf) error "Set --gtf (needed for gene names and counting)"
    if (!(params.strandedness in ['auto', 'unstranded', 'forward', 'reverse']))
        error "--strandedness must be auto, unstranded, forward or reverse"

    def sheet = file(params.samplesheet, checkIfExists: true)
    def base  = sheet.parent

    ch_reads = Channel.fromPath(sheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.condition || !row.fastq_1)
                error "Samplesheet row missing sample, condition or fastq_1: ${row}"
            def meta  = [id: row.sample, condition: row.condition]
            def reads = [resolvePath(row.fastq_1, base)]
            if (row.fastq_2) reads << resolvePath(row.fastq_2, base)
            tuple(meta, reads)
        }

    gtf = file(params.gtf, checkIfExists: true)
    index = params.star_index
        ? Channel.value(file(params.star_index, checkIfExists: true))
        : STAR_INDEX(file(params.genome, checkIfExists: true), gtf).index

    FASTQC(ch_reads)
    STAR_ALIGN(ch_reads, index)

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
