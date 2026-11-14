// Align reads with STAR and count reads per gene in the same pass (--quantMode GeneCounts).
// Several runs of one sample are aligned together (comma-separated lists for STAR);
// paired-end samples pass read 1 and read 2 lists separately.
process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::star=2.7.11b'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'

    input:
    tuple val(meta), path(r1), path(r2)
    path index

    output:
    tuple val(meta), path("${meta.id}.ReadsPerGene.out.tab"),          emit: counts
    tuple val(meta), path("${meta.id}.Log.final.out"),                 emit: log
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"), emit: bam
    path 'versions.yml', emit: versions

    script:
    def r1s = [r1].flatten().collect { it.name }
    def r2s = [r2].flatten().collect { it.name }
    def gz  = r1s[0].endsWith('.gz')
    // The STAR biocontainer ships several CPU-specific builds and picks one from /proc/cpuinfo.
    // On mixed-CPU clusters that choice can crash (SIGILL); --star_simd forces one build, e.g. sse4.1.
    def star = params.star_simd ? "STAR-${params.star_simd}" : 'STAR'
    // STAR streams .gz input through --readFilesCommand on Linux (HPC, CI). STAR builds for
    // macOS cannot spawn that command, so on macOS the reads are decompressed first.
    """
    R1="${r1s.join(' ')}"
    R2="${r2s.join(' ')}"
    CMD=""
    if [ "${gz}" = "true" ]; then
        if [ "\$(uname)" = "Darwin" ]; then
            unz() { out=""; for f in "\$@"; do gunzip -c "\$f" > "\${f%.gz}.unzipped"; out="\$out \${f%.gz}.unzipped"; done; echo \$out; }
            R1=\$(unz \$R1)
            if [ -n "\$R2" ]; then R2=\$(unz \$R2); fi
        else
            CMD="--readFilesCommand gunzip -c"
        fi
    fi
    J1=\$(echo \$R1 | tr ' ' ',')
    J2=\$(echo \$R2 | tr ' ' ',')
    ${star} --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn \$J1 \$J2 \\
        \$CMD \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode GeneCounts \\
        --outFileNamePrefix ${meta.id}.
    rm -f *.unzipped
    printf '"%s":\\n    star: "%s"\\n' "${task.process}" "\$(${star} --version)" > versions.yml
    """
}
