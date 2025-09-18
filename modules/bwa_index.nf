process BWA_INDEX {
    label 'process_high'
    
    publishDir "${params.outdir}/genome", mode: 'copy'
    
    input:
    path fasta
    
    output:
    tuple val("${fasta.baseName}"), path("bwa_index"), emit: index
    path "versions.yml", emit: versions
    
    script:
    """
    mkdir bwa_index
    bwa index -p bwa_index/${fasta.baseName} $fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}