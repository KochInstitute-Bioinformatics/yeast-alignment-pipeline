include { FASTQC } from '../modules/fastqc'
include { CONCATENATE_GENOME } from '../modules/concatenate_genome'
include { BWA_INDEX } from '../modules/bwa_index'
include { BWA_ALIGN } from '../modules/bwa_align'
include { SAMTOOLS_SORT } from '../modules/samtools_sort'
include { SAMTOOLS_INDEX } from '../modules/samtools_index'
include { SAMTOOLS_FAIDX } from '../modules/samtools_faidx'
include { MULTIQC } from '../modules/multiqc'

workflow YEAST_ALIGNMENT {
    
    take:
    input_csv
    reference_fasta
    
    main:
    // Parse input samplesheet
    reads_ch = Channel
        .fromPath(input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.single_end = row.fastq_2 ? false : true
            meta.vector = row.vector ?: null
            
            def reads = meta.single_end ? 
                [file(row.fastq_1, checkIfExists: true)] :
                [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
            
            return [meta, reads]
        }
    
    // Debug: Print number of samples
    reads_ch.count().view { "Total samples found: $it" }
    
    // Create reference genome channel
    reference_ch = Channel.fromPath(reference_fasta, checkIfExists: true)
    
    // Group samples by their vector file to create unique reference combinations
    vector_groups = reads_ch
        .map { meta, reads -> 
            def vector_key = meta.vector ? file(meta.vector).baseName : "WT_only"
            def vector_file = meta.vector ? file(meta.vector, checkIfExists: true) : null
            return [vector_key, vector_file, meta, reads]
        }
        .groupTuple(by: [0, 1])
        .map { vector_key, vector_file_list, metas, reads_list ->
            // Handle the vector file - take first non-null entry
            def vector_file = vector_file_list.find { it != null }
            return [vector_key, vector_file, metas, reads_list]
        }
    
    // Create unique reference combinations
    reference_combinations = vector_groups
        .map { vector_key, vector_file, metas, reads_list ->
            return [vector_key, vector_file]
        }
        .unique()
    
    // Create combined references for each unique combination
    CONCATENATE_GENOME(
        reference_combinations.combine(reference_ch)
            .map { vector_key, vector_file, ref_fasta -> 
                def vector_files = vector_file ? [vector_file] : []
                return [vector_key, ref_fasta, vector_files]
            }
    )
    
    // Create FASTA index for each reference
    SAMTOOLS_FAIDX(CONCATENATE_GENOME.out.fasta)
    
    // Index each unique reference with BWA
    BWA_INDEX(CONCATENATE_GENOME.out.fasta)
    
    // Run FastQC on raw reads
    if (!params.skip_fastqc) {
        FASTQC(reads_ch)
        fastqc_ch = FASTQC.out.zip.map { meta, files -> files }.flatten()
    } else {
        fastqc_ch = Channel.empty()
    }
    
    // Prepare alignment input by matching samples to their appropriate reference
    alignment_input = vector_groups
        .transpose()
        .map { vector_key, vector_file, meta, reads ->
            return [vector_key, meta, reads]
        }
        .combine(BWA_INDEX.out.index, by: 0)
        .map { vector_key, meta, reads, index_name, index_path ->
            return [meta, reads, index_name, index_path]
        }
    
    // Align reads to appropriate reference
    BWA_ALIGN(alignment_input)
    
    // Sort BAM files
    SAMTOOLS_SORT(BWA_ALIGN.out.bam)
    
    // Index sorted BAM files
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    
    // Collect all QC files for MultiQC
    if (!params.skip_multiqc) {
        multiqc_files = Channel.empty()
            .mix(fastqc_ch)
            .mix(BWA_ALIGN.out.log.map { meta, files -> files }.flatten())
            .mix(SAMTOOLS_SORT.out.stats.map { meta, files -> files }.flatten())
            .collect()
        
        MULTIQC(multiqc_files)
    }
    
    emit:
    bam = SAMTOOLS_SORT.out.bam
    bai = SAMTOOLS_INDEX.out.bai
    combined_references = CONCATENATE_GENOME.out.fasta
    reference_indices = SAMTOOLS_FAIDX.out.fai
}