# Automated Reference-Based Assembly Pipeline with Host Read Removal

## Overview

This pipeline automates the process of reference-based assembly for genomic sequencing data, including the removal of host reads using **Hostile**. It is designed to streamline the workflow from raw paired-end sequencing reads to consensus genome generation. The script integrates multiple tools for quality control, read trimming, host read filtering, alignment, variant calling, and consensus sequence construction.

### Key Features:
- **Host Read Removal**: Filters host reads using **Hostile** to improve downstream analysis.
- **Reference-Based Assembly**: Aligns filtered reads to a reference genome for variant calling and consensus sequence generation.
- **Comprehensive Workflow**: Combines quality control, alignment, and variant analysis into a single automated process.
- **High-Throughput**: Processes multiple samples in parallel using a sample list.

---

## Repository

This repository is hosted on GitHub:  
[Automated-Reference-Based-Assembly-Pipeline-with-Host-Read-Removal](https://github.com/gmboowa/Automated-Reference-Based-Assembly-Pipeline-with-Host-Read-Removal)

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
- **Reference Genome**: A FASTA file of the target genome (e.g., viral genome).
- **Host Reference Genome**: A FASTA file of the host genome (e.g., human genome for filtering host reads).
- **Paired-End Sequencing Reads**: Raw FASTQ files for each sample named as:
  - `sampleID_1.clean_1.fastq.gz`
  - `sampleID_2.clean_2.fastq.gz`
- **Sample List**: A text file (`samples.txt`) with one sample ID per line, e.g.:
  ```
  sample1
  sample2
  sample3

Usage
1. Clone the Repository
To get started, you’ll need to clone this repository from GitHub. Run the following commands in your terminal:


git clone https://github.com/gmboowa/Automated-Reference-Based-Assembly-Pipeline-with-Host-Read-Removal.git
cd Automated-Reference-Based-Assembly-Pipeline-with-Host-Read-Removal

2. Make the Script Executable
To ensure the script can run properly, make it executable by using:


chmod +x reference_assembly.sh
3. Prepare Input Files
Before running the pipeline, make sure you have all the required input files organized in your working directory:

Reference Genome: A FASTA file of the target genome (e.g., reference.fasta).
Host Reference Genome: A FASTA file of the host genome (e.g., host_reference.fasta).
Paired-End Reads: Raw FASTQ files for each sample, named as:
sampleID_1.clean_1.fastq.gz
sampleID_2.clean_2.fastq.gz
Sample List: A text file (samples.txt) with one sample ID per line:

sample1
sample2
sample3

4. Run the Pipeline
Once everything is ready, you can run the pipeline to process your samples by simply executing:


./reference_assembly.sh

Workflow Steps
This pipeline follows a structured workflow:

Indexing
Indexes the reference genome (target) and host genome to prepare for alignment and filtering.

Quality Control
Runs FastQC on raw sequencing reads to evaluate their quality.

Read Trimming
Trims adapters and low-quality regions using fastp.

Host Read Removal
Filters out host reads from the dataset using Hostile.

Alignment
Aligns the host-removed reads to the reference genome using either BWA or Minimap2.

Variant Calling
Identifies variants in the aligned reads using BCFtools.

Consensus Genome Construction
Creates a consensus genome based on the identified variants using BCFtools consensus.

Coverage Analysis (Optional)
Calculates genome coverage using Samtools.

Example Outputs
After processing each sample, the pipeline generates the following outputs:

Quality-Control Reports:
sample1_1_fastqc.html, sample1_2_fastqc.html
Trimmed Reads:
sample1_trimmed_1.fastq.gz, sample1_trimmed_2.fastq.gz
Host-Filtered Reads:
sample1_host_removed_r1.fastq.gz, sample1_host_removed_r2.fastq.gz
Alignment Files:
sample1_sorted.bam, sample1_sorted.bam.bai
Variant Files:
sample1_variants.vcf.gz, sample1_variants.vcf.gz.tbi
Consensus Genome:
sample1_consensus.fasta
Coverage File:
sample1_coverage.txt
Troubleshooting
Missing Dependencies
If any required tool is missing, you can easily install it using Conda:


conda install -c bioconda bwa minimap2 fastqc fastp samtools bcftools
Script Fails Midway
Check the log file generated for each sample for details on any errors.
Ensure all input files are correctly named and located in the working directory.
License
This pipeline is open-source and distributed under the MIT License. For details, see the LICENSE file included in this repository.

Acknowledgments
This pipeline integrates tools developed by the open-source bioinformatics community, including BWA, Minimap2, Hostile, FastQC, fastp, Samtools, and BCFtools. Special thanks to the developers for their contributions!



This is the complete and polished **README.md** file formatted for GitHub. You can copy and paste this directly into your repository. Let me know if you’d like additional changes!













