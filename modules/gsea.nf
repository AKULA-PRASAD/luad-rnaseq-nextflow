// Pre-ranked gene set enrichment (fgsea) on the DESeq2 Wald statistic.
// Default gene sets: MSigDB Hallmark (human gene symbols, CC BY 4.0), bundled in assets/.
process GSEA {
    label 'process_low'
    conda 'bioconda::bioconductor-fgsea=1.32.2 conda-forge::r-base=4.4'
    container 'quay.io/biocontainers/bioconductor-fgsea:1.32.2--r44he5774e6_0'

    input:
    path results
    path gmt

    output:
    path 'gsea_results.tsv', emit: results
    path 'gsea_summary.txt', emit: summary
    path 'gsea_nes.pdf',     emit: plot
    path 'versions.yml',     emit: versions

    script:
    """
    gsea.R --results ${results} --gmt ${gmt} --padj ${params.padj} \\
        --min_size ${params.gsea_min_size} --max_size ${params.gsea_max_size} --seed 1
    printf '"%s":\\n    r: "%s"\\n    fgsea: "%s"\\n' "${task.process}" \\
        "\$(Rscript -e 'cat(as.character(getRversion()))')" \\
        "\$(Rscript -e 'cat(as.character(packageVersion("fgsea")))')" > versions.yml
    """
}
