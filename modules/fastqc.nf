process FASTQC {
    label 'process_medium'
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastqc $args --threads $task.cpus $reads
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | tail -1 | sed 's/FastQC v//')
    END_VERSIONS
    """
}