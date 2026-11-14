// Build the STAR genome index once; reused by every alignment task.
// sjdbOverhang = read length - 1 (77 for the 2 x 78 bp reads of GSE37764; 26 for the 27 bp test reads).
process STAR_INDEX {
    tag "${fasta.simpleName}"
    label 'process_high'
    conda 'bioconda::star=2.7.11b'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'

    input:
    path fasta
    path gtf

    output:
    path 'star_index', emit: index
    path 'versions.yml', emit: versions

    script:
    def fa_gz  = fasta.name.endsWith('.gz')
    def gtf_gz = gtf.name.endsWith('.gz')
    def fa     = fa_gz  ? fasta.baseName : fasta.name
    def ann    = gtf_gz ? gtf.baseName   : gtf.name
    def ram = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes()}" : ''
    """
    ${fa_gz  ? "gunzip -c ${fasta} > ${fa}" : ''}
    ${gtf_gz ? "gunzip -c ${gtf} > ${ann}"  : ''}
    mkdir star_index
    STAR --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fa} \\
        --sjdbGTFfile ${ann} \\
        --sjdbOverhang ${params.sjdb_overhang} \\
        --genomeSAindexNbases ${params.sa_index_nbases} \\
        ${ram}
    printf '"%s":\\n    star: "%s"\\n' "${task.process}" "\$(STAR --version)" > versions.yml
    """
}
