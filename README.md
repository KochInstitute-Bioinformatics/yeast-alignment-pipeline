# Yeast DNA Alignment Pipeline with Vector Detection

A Nextflow pipeline for aligning DNA sequencing reads to yeast genomes with optional vector sequence integration.

## Features

- **FastQC**: Quality control of raw sequencing reads
- **Vector Integration**: Concatenates vector sequences to reference genome
- **BWA Alignment**: Aligns reads to combined reference (genome + vectors)
- **BAM Processing**: Sorts and indexes alignment files
- **MultiQC**: Comprehensive quality control reporting
- **Singularity Support**: Containerized execution for reproducibility

## Quick Start

### Prerequisites

- Nextflow (>=23.04.0)
- Singularity
- Linux environment

### Installation

```bash
git clone <your-repo-url>
cd yeast_alignment_pipeline
```

### Usage

1. **Prepare your samplesheet** (`samplesheet.csv`):
```bash
sample,fastq_1,fastq_2,vector
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,/path/to/vector1.fasta
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,
sample3,/path/to/sample3_R1.fastq.gz,,/path/to/vector2.fasta
```

2. **Run the pipeline**:
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --fasta /path/to/yeast_genome.fasta \
    --outdir results
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to samplesheet CSV | Required |
| `--fasta` | Path to reference genome FASTA | Required |
| `--outdir` | Output directory | `results` |
| `--skip_fastqc` | Skip FastQC step | `false` |
| `--skip_multiqc` | Skip MultiQC step | `false` |

### Output Structure

```bash
results/
├── alignments/           # Sorted BAM files and indices
├── fastqc/              # FastQC reports
├── genome/              # Combined reference genome
├── multiqc/             # MultiQC report
└── pipeline_info/       # Execution reports
```

## Pipeline Overview

1. **Input Parsing**: Reads samplesheet and validates file paths
2. **Vector Integration**: Concatenates unique vector sequences to genome
3. **Reference Indexing**: Creates BWA index for combined reference
4. **Quality Control**: Runs FastQC on raw reads
5. **Alignment**: Aligns reads using BWA-MEM
6. **Post-processing**: Sorts and indexes BAM files
7. **Reporting**: Generates comprehensive MultiQC report

## Vector Handling

- Vector sequences are concatenated to the reference genome with modified headers (`vector_<name>_<sequence>`)
- All samples are aligned to the same combined reference containing all unique vectors
- Reads aligning to `vector_*` sequences indicate vector presence/contamination
- Samples without vectors simply won't have reads aligning to vector sequences

## Configuration

The pipeline uses Singularity containers by default. Modify `nextflow.config` to adjust:
- Resource requirements
- Container settings
- Process-specific parameters

## Troubleshooting

### Common Issues

1. **Container pull failures**: Increase `singularity.pullTimeout` in config
2. **Memory errors**: Adjust memory settings in process labels
3. **File not found**: Ensure all paths in samplesheet are absolute paths

### Resume Failed Runs

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --fasta genome.fasta \
    --outdir results \
    -resume
```

## Citation

If you use this pipeline, please cite:
- [Nextflow](https://www.nextflow.io/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [Samtools](http://www.htslib.org/)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

- Author: Charlie W
- Institution: [Your Institution]
- Email: [Your Email]
