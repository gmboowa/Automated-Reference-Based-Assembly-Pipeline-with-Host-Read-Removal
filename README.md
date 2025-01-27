# Automated Reference-Based Assembly Pipeline with Host Read Removal

## Overview

This pipeline automates the process of reference-based assembly for genomic sequencing data, including the removal of host reads using **Hostile**. It is designed to streamline the workflow from raw paired-end sequencing reads to consensus genome generation. The script integrates multiple tools for quality control, read trimming, host read filtering, alignment, variant calling, and consensus sequence construction.

### Key Features:
- **Host Read Removal**: Filters host reads using **Hostile** to improve downstream analysis.
- **Reference-Based Assembly**: Aligns filtered reads to a reference genome for variant calling and consensus sequence generation.
- **Comprehensive Workflow**: Combines quality control, alignment, and variant analysis into a single automated process.
- **High-Throughput**: Processes multiple samples in parallel using a sample list.

---

## Alternative Names

This pipeline can also be referred to as:
- **Comprehensive Pipeline for Reference-Based Assembly with Hostile Integration**
- **Host-Read Removal and Viral Reference Assembly Workflow**
- **Reference-Based Assembly of Viral Genomes with Host Read Filtering**
- **High-Throughput Assembly Pipeline: Host Filtering + Variant Analysis**

---

## Prerequisites

### Required Tools
Ensure the following tools are installed and accessible in your `$PATH`:
- [**BWA**](http://bio-bwa.sourceforge.net/)
- [**Minimap2**](https://github.com/lh3/minimap2)
- [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [**fastp**](https://github.com/OpenGene/fastp)
- [**Hostile**](https://github.com/bede/hostile)
- [**Samtools**](http://www.htslib.org/)
- [**BCFtools**](http://www.htslib.org/)

### Required Files
- **Reference Genome**: FASTA file of the target genome (e.g., viral genome).
- **Host Reference Genome**: FASTA file of the host genome (e.g., human genome for filtering host reads).
- **Paired-End Sequencing Reads**: Raw FASTQ files named as:
  - `sampleID_1.clean_1.fastq.gz`
  - `sampleID_2.clean_2.fastq.gz`
- **Sample List**: A text file (`samples.txt`) containing one sample ID per line:
  ```text
  sample1
  sample2
  sample3
