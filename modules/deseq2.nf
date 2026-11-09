// Tumour vs normal differential expression (DESeq2: median-of-ratios, Wald test, BH).
process DESEQ2 {
    label 'process_medium'
    conda 'bioconda::bioconductor-deseq2=1.46.0 conda-forge::r-base=4.4'
    container 'quay.io/biocontainers/bioconductor-deseq2:1.46.0--r44he5774e6_1'

    input:
    path counts
    path samplesheet

    output:
    path 'deseq2_results.csv',  emit: results
    path 'significant_genes.csv', emit: significant
    path 'summary.txt',          emit: summary
    path '*.pdf',                emit: plots

    script:
    """
    deseq2.R --counts ${counts} --samplesheet ${samplesheet} \\
        --reference ${params.reference_level} --padj ${params.padj}
    """
}
