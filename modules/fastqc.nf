// Read-level quality control for each sample.
process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::fastqc=0.12.1'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.html'), emit: html
    tuple val(meta), path('*.zip'),  emit: zip

    // FastQC is a Java program. Each JVM writes /tmp/hsperfdata_<user>/<pid>; parallel tasks in
    // separate containers can share a pid and the host /tmp, collide on that file and crash with
    // SIGBUS. -XX:-UsePerfData turns the file off; --dir keeps FastQC's temp files in the work dir.
    script:
    """
    export _JAVA_OPTIONS="-XX:-UsePerfData"
    fastqc --threads ${task.cpus} --quiet --dir . ${reads}
    """
}
