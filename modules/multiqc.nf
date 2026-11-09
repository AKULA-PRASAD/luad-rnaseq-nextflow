// One HTML report combining FastQC and STAR alignment statistics for all samples.
process MULTIQC {
    label 'process_low'
    conda 'bioconda::multiqc=1.33'
    container 'quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0'

    input:
    path files

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data

    script:
    """
    multiqc --force --filename multiqc_report.html .
    """
}
