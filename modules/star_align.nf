// Align reads with STAR and count reads per gene in the same pass (--quantMode GeneCounts).
process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::star=2.7.11b'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}.ReadsPerGene.out.tab"),      emit: counts
    tuple val(meta), path("${meta.id}.Log.final.out"),             emit: log
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"), emit: bam

    script:
    def files = [reads].flatten()
    def gz    = files[0].name.endsWith('.gz')
    // The STAR biocontainer ships several CPU-specific builds and picks one from /proc/cpuinfo.
    // On mixed-CPU clusters that choice can crash (SIGILL); --star_simd forces one build, e.g. sse4.1.
    def star  = params.star_simd ? "STAR-${params.star_simd}" : 'STAR'
    // STAR streams .gz input through --readFilesCommand on Linux (HPC, CI). STAR builds for
    // macOS cannot spawn that command, so on macOS the reads are decompressed first.
    """
    READS="${files.join(' ')}"
    CMD=""
    if [ "${gz}" = "true" ]; then
        if [ "\$(uname)" = "Darwin" ]; then
            READS=""
            for f in ${files.join(' ')}; do gunzip -c "\$f" > "\${f%.gz}.unzipped"; READS="\$READS \${f%.gz}.unzipped"; done
        else
            CMD="--readFilesCommand gunzip -c"
        fi
    fi
    ${star} --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn \$READS \\
        \$CMD \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode GeneCounts \\
        --outFileNamePrefix ${meta.id}.
    rm -f *.unzipped
    """
}
