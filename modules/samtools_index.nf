process SAMTOOLS_INDEX {
    conda "bioconda::samtools=1.21"
    label 'process_low'
    
    publishDir "${params.outdir}/alignments", mode: 'copy'
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    """
    samtools index $args $bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}