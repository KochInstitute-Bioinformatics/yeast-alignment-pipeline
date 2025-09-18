process MULTIQC {
    label 'process_low'
    
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path multiqc_files
    
    output:
    path "*.html", emit: report
    path "*_data", emit: data, optional: true
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    """
    multiqc $args .
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(echo \$(multiqc --version 2>&1) | sed 's/^.*multiqc, version //; s/ .*\$//')
    END_VERSIONS
    """
}