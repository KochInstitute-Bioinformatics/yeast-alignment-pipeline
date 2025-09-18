process CONCATENATE_GENOME {
    label 'process_low'
    
    publishDir "${params.outdir}/genome", mode: 'copy'
    
    input:
    path genome_fasta
    path vector_files
    
    output:
    path "combined_reference.fasta", emit: fasta
    path "reference_info.txt", emit: info
    path "versions.yml", emit: versions
    
    script:
    def vectors_present = vector_files.size() > 0
    """
    # Start with the original genome
    cp $genome_fasta combined_reference.fasta
    
    # Create info file
    echo "Reference genome: $genome_fasta" > reference_info.txt
    echo "Original genome sequences:" >> reference_info.txt
    grep "^>" $genome_fasta | wc -l >> reference_info.txt
    
    if [ "$vectors_present" = "true" ]; then
        echo "Vector files added:" >> reference_info.txt
        
        # Add each vector file to the combined reference
        for vector_file in $vector_files; do
            echo "Adding vector: \$vector_file" >> reference_info.txt
            
            # Modify headers to include vector identifier
            vector_name=\$(basename \$vector_file .fasta)
            vector_name=\$(basename \$vector_name .fa)
            
            # Add vector sequences with modified headers
            awk -v vname="\$vector_name" '
                /^>/ { 
                    gsub(/^>/, ">vector_" vname "_")
                    print
                } 
                !/^>/ { print }
            ' \$vector_file >> combined_reference.fasta
            
            echo "  - \$vector_file (as vector_\${vector_name}_*)" >> reference_info.txt
        done
    else
        echo "No vector files provided - using original genome only" >> reference_info.txt
    fi
    
    # Final statistics
    echo "" >> reference_info.txt
    echo "Final combined reference statistics:" >> reference_info.txt
    echo "Total sequences: \$(grep '^>' combined_reference.fasta | wc -l)" >> reference_info.txt
    echo "Total length: \$(grep -v '^>' combined_reference.fasta | tr -d '\\n' | wc -c)" >> reference_info.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | sed 's/^.*awk //; s/,.*\$//')
    END_VERSIONS
    """
}