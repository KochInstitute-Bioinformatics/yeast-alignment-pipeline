#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YEAST_ALIGNMENT } from './workflows/yeast_alignment'

workflow {
    // Check required parameters
    if (!params.input) {
        error "Please provide an input samplesheet with --input"
    }
    if (!params.fasta) {
        error "Please provide a reference genome with --fasta"
    }
    
    YEAST_ALIGNMENT(params.input, params.fasta)
}