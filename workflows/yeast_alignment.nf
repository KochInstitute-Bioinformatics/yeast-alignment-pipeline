include { FASTQC } from '../modules/fastqc'
include { CONCATENATE_GENOME } from '../modules/concatenate_genome'
include { BWA_INDEX } from '../modules/bwa_index'
include { BWA_ALIGN } from '../modules/bwa_align'
include { SAMTOOLS_SORT } from '../modules/samtools_sort'
include { SAMTOOLS_INDEX } from '../modules/samtools_index'
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
    
    // Collect all unique vector files
    vector_files_ch = reads_ch
        .filter { meta, reads -> meta.vector != null && meta.vector != "" }
        .map { meta, reads -> file(meta.vector, checkIfExists: true) }
        .unique()
        .collect()
        .ifEmpty([])
    
    // Concatenate genome with vectors (if any vectors exist)
    CONCATENATE_GENOME(reference_ch, vector_files_ch)
    
    // Run FastQC on raw reads (this should process all samples)
    if (!params.skip_fastqc) {
        FASTQC(reads_ch)
        fastqc_ch = FASTQC.out.zip.map { meta, files -> files }.flatten()
    } else {
        fastqc_ch = Channel.empty()
    }
    
    // Index the combined reference genome
    BWA_INDEX(CONCATENATE_GENOME.out.fasta)
    
    // Align reads to combined reference - FIX: Use combine instead of passing tuple
    BWA_ALIGN(reads_ch.combine(BWA_INDEX.out.index))
    
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
    combined_reference = CONCATENATE_GENOME.out.fasta
}