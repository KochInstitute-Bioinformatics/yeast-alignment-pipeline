process BWA_ALIGN {
    label 'process_high'
    
    input:
    tuple val(meta), path(reads), val(index_name), path(index)
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA"
    
    if (meta.single_end) {
        """
        bwa mem \\
            $args \\
            -t $task.cpus \\
            -R '$read_group' \\
            ${index}/${index_name} \\
            $reads \\
            2> ${prefix}.bwa.log \\
            | samtools view -@ $task.cpus -b -h -o ${prefix}.bam -
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        bwa mem \\
            $args \\
            -t $task.cpus \\
            -R '$read_group' \\
            ${index}/${index_name} \\
            ${reads[0]} ${reads[1]} \\
            2> ${prefix}.bwa.log \\
            | samtools view -@ $task.cpus -b -h -o ${prefix}.bam -
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}