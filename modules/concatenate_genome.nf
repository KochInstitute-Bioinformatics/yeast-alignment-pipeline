process CONCATENATE_GENOME {
    label 'process_low'
    publishDir "${params.outdir}/genome", mode: 'copy'
    
    input:
    tuple val(vector_key), path(genome_fasta), path(vector_files)
    
    output:
    tuple val(vector_key), path("${vector_key}.fasta"), emit: fasta
    path "${vector_key}_info.txt", emit: info
    path "versions.yml", emit: versions
    
    script:
    def vectors_present = vector_files.size() > 0
    def output_name = "${vector_key}.fasta"
    """
    # Start with the original genome
    cp $genome_fasta $output_name
    
    # Create info file
    echo "Reference: $vector_key" > ${vector_key}_info.txt
    echo "Base genome: $genome_fasta" >> ${vector_key}_info.txt
    echo "Original genome sequences: \$(grep '^>' $genome_fasta | wc -l)" >> ${vector_key}_info.txt
    
    if [ "$vectors_present" = "true" ]; then
        echo "Vector files added:" >> ${vector_key}_info.txt
        
        # Add each vector file to the combined reference
        for vector_file in $vector_files; do
            echo "Adding vector: \$vector_file" >> ${vector_key}_info.txt
            
            # Add vector sequences with original FASTA headers (just the defline value)
            cat \$vector_file >> $output_name
            
            # Extract the header names for logging
            vector_headers=\$(grep '^>' \$vector_file | sed 's/^>//' | tr '\\n' ',' | sed 's/,\$//')
            echo "  - \$vector_file (headers: \$vector_headers)" >> ${vector_key}_info.txt
        done
    else
        echo "No vector files provided - using WT genome only" >> ${vector_key}_info.txt
    fi
    
    # Final statistics
    echo "" >> ${vector_key}_info.txt
    echo "Final combined reference statistics:" >> ${vector_key}_info.txt
    echo "Total sequences: \$(grep '^>' $output_name | wc -l)" >> ${vector_key}_info.txt
    echo "Total length: \$(grep -v '^>' $output_name | tr -d '\\n' | wc -c)" >> ${vector_key}_info.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | sed 's/^.*awk //; s/,.*\$//')
    END_VERSIONS
    """
}