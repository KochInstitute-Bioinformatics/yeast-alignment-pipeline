process BWA_INDEX {
    label 'process_high'
    
    publishDir "${params.outdir}/genome", mode: 'copy'
    
    input:
    tuple val(vector_key), path(fasta)
    
    output:
    tuple val(vector_key), val("${fasta.baseName}"), path("bwa_index_${vector_key}"), emit: index
    path "versions.yml", emit: versions
    
    script:
    """
    mkdir bwa_index_${vector_key}
    bwa index -p bwa_index_${vector_key}/${fasta.baseName} $fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}