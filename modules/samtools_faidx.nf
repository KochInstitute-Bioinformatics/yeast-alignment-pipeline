process SAMTOOLS_FAIDX {
    label 'process_low'
    
    publishDir "${params.outdir}/genome", mode: 'copy'
    
    input:
    tuple val(vector_key), path(fasta)
    
    output:
    tuple val(vector_key), path("*.fai"), emit: fai
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    """
    samtools faidx $args $fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}