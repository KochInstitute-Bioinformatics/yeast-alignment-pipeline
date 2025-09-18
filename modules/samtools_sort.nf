process SAMTOOLS_SORT {
    conda "bioconda::samtools=1.21"
    label 'process_medium'
    
    publishDir "${params.outdir}/alignments", mode: 'copy'
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    tuple val(meta), path("*.stats"), emit: stats
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.sorted.bam \\
        $bam
    
    samtools stats ${prefix}.sorted.bam > ${prefix}.stats
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}