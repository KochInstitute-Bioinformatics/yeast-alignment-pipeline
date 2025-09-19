process BAMCOVERAGE_BIGWIG {
    container "quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1"
    label 'process_medium'
    publishDir "${params.outdir}/alignments", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: '--normalizeUsing CPM --binSize 10'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}.bw \\
        --outFileFormat bigwig \\
        --numberOfProcessors $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
    END_VERSIONS
    """
}