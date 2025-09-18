#!/bin/bash

# Example run script for yeast alignment pipeline

# Set paths (modify these for your data)
SAMPLESHEET="examples/samplesheet.csv"
GENOME="/path/to/yeast_genome.fasta"
OUTDIR="results"

# Run pipeline
nextflow run main.nf \
    --input $SAMPLESHEET \
    --fasta $GENOME \
    --outdir $OUTDIR \
    -resume

echo "Pipeline completed! Check results in: $OUTDIR"
