// Optional adapter and quality trimming (--trim). Off by default: STAR soft-clips poor read ends,
// so trimming is mainly useful for libraries with heavy adapter read-through.
// Runs per sequencing run; the trimmed reads go to STAR, the JSON report to MultiQC.
process FASTP {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::fastp=1.3.6'
    container 'quay.io/biocontainers/fastp:1.3.6--h43da1c4_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.fastp.json'),    emit: json
    tuple val(meta), path('*.fastp.html'),    emit: html
    path 'versions.yml', emit: versions

    script:
    def files = [reads].flatten()
    def io = files.size() == 2
        ? "--in1 ${files[0]} --in2 ${files[1]} --out1 ${meta.id}_1.trim.fastq.gz --out2 ${meta.id}_2.trim.fastq.gz --detect_adapter_for_pe"
        : "--in1 ${files[0]} --out1 ${meta.id}.trim.fastq.gz"
    """
    fastp ${io} \\
        --thread ${task.cpus} \\
        --qualified_quality_phred ${params.trim_min_quality} \\
        --length_required ${params.trim_min_length} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    printf '"%s":\\n    fastp: "%s"\\n' "${task.process}" "\$(fastp --version 2>&1 | sed 's/fastp //')" > versions.yml
    """
}
