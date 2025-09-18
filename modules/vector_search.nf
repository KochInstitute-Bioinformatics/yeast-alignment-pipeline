process VECTOR_SEARCH {
    label 'process_medium'
    
    publishDir "${params.outdir}/vector_analysis", mode: 'copy'
    
    input:
    tuple val(meta), path(bam), path(vector_fasta)
    tuple val(index_name), path(index)
    
    output:
    tuple val(meta), path("*.vector_hits.txt"), emit: results
    tuple val(meta), path("*.vector_reads.fastq"), emit: reads, optional: true
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_identity = params.vector_min_identity ?: 90
    def min_length = params.vector_min_length ?: 50
    """
    # Extract unmapped reads (potential vector insertions)
    samtools view -f 4 -h $bam | samtools fastq - > ${prefix}.unmapped.fastq
    
    # Create BLAST database from vector sequence
    makeblastdb -in $vector_fasta -dbtype nucl -out vector_db
    
    # Search for vector sequences in unmapped reads
    blastn \\
        -query ${prefix}.unmapped.fastq \\
        -db vector_db \\
        -out ${prefix}.vector_hits.txt \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
        -perc_identity $min_identity \\
        -qcov_hsp_perc 80 \\
        -num_threads $task.cpus
    
    # Extract reads with significant vector hits
    if [ -s ${prefix}.vector_hits.txt ]; then
        cut -f1 ${prefix}.vector_hits.txt | sort -u > vector_read_ids.txt
        
        # Check if seqtk is available, otherwise use awk
        if command -v seqtk &> /dev/null; then
            seqtk subseq ${prefix}.unmapped.fastq vector_read_ids.txt > ${prefix}.vector_reads.fastq
        else
            # Alternative method using awk if seqtk is not available
            awk 'NR==FNR{ids[\$1]; next} /^@/{f=(\$1 in ids)} f' vector_read_ids.txt ${prefix}.unmapped.fastq > ${prefix}.vector_reads.fastq
        fi
        
        echo "Found \$(wc -l < vector_read_ids.txt) reads with vector sequences" >> ${prefix}.vector_hits.txt
    else
        echo "No vector sequences found in unmapped reads" > ${prefix}.vector_hits.txt
    fi
    
    # Clean up
    rm -f ${prefix}.unmapped.fastq vector_db.* vector_read_ids.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}