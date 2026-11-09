// Merge per-sample STAR gene counts into one matrix, choosing the strand column
// automatically (or as set by --strandedness) and adding gene names from the GTF.
process MERGE_COUNTS {
    label 'process_low'
    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    path tabs
    path gtf
    val strandedness

    output:
    path 'counts_matrix.csv',    emit: matrix
    path 'strandedness.tsv',     emit: strandedness

    script:
    """
    merge_counts.py --gtf ${gtf} --strandedness ${strandedness} \\
        --out counts_matrix.csv --report strandedness.tsv ${tabs}
    """
}
