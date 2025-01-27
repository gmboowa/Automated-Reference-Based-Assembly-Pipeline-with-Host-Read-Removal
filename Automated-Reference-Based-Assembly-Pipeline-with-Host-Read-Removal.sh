#!/bin/bash

# Exit script if any command fails
set -e

# Input variables
SAMPLES_LIST="samples.txt"    # File containing sample IDs (one per line)
REFERENCE="reference.fasta"   # Reference genome file for the target (e.g., virus)
HOST_REFERENCE="host_reference.fasta"  # Host reference genome file for Hostile
THREADS=4                     # Number of threads to use

# Ensure required tools are installed
for tool in bwa minimap2 fastqc fastp samtools bcftools hostile; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool is not installed. Please install it before running the script."
        exit 1
    fi
done

# Check if the reference genome exists
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference genome file '$REFERENCE' not found."
    exit 1
fi

# Check if the host reference genome exists
if [ ! -f "$HOST_REFERENCE" ]; then
    echo "Error: Host reference genome file '$HOST_REFERENCE' not found."
    exit 1
fi

# Check if the samples list exists
if [ ! -f "$SAMPLES_LIST" ]; then
    echo "Error: Samples list file '$SAMPLES_LIST' not found."
    exit 1
fi

# Step 1: Index the reference genomes
echo "Indexing the target reference genome..."
bwa index $REFERENCE
minimap2 -d reference.mmi $REFERENCE

echo "Indexing the host reference genome for Hostile..."
bwa index $HOST_REFERENCE

# Loop through each sample
while read -r SAMPLE_ID; do
    echo "Processing sample: $SAMPLE_ID"

    # Define file paths for paired-end reads
    R1="${SAMPLE_ID}_1.clean_1.fastq.gz"
    R2="${SAMPLE_ID}_2.clean_2.fastq.gz"

    # Check if raw reads exist
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "Error: Raw read files for sample '$SAMPLE_ID' not found."
        echo "Missing files: $R1 or $R2"
        continue  # Skip to the next sample
    fi

    # Step 2: Quality control (Optional)
    echo "Running FastQC for $SAMPLE_ID..."
    fastqc $R1 $R2

    # Step 3: Trim adapters and low-quality bases
    echo "Trimming reads for $SAMPLE_ID..."
    fastp -i $R1 -I $R2 -o ${SAMPLE_ID}_trimmed_1.fastq.gz -O ${SAMPLE_ID}_trimmed_2.fastq.gz -w $THREADS

    # Update trimmed file paths
    TRIMMED_R1="${SAMPLE_ID}_trimmed_1.fastq.gz"
    TRIMMED_R2="${SAMPLE_ID}_trimmed_2.fastq.gz"

    # Step 4: Remove host reads using Hostile
    echo "Removing host reads for $SAMPLE_ID..."
    hostile --host-index $HOST_REFERENCE \
            --threads $THREADS \
            --r1 $TRIMMED_R1 \
            --r2 $TRIMMED_R2 \
            --output-prefix ${SAMPLE_ID}_host_removed

    # Update file paths after Hostile filtering
    HOST_REMOVED_R1="${SAMPLE_ID}_host_removed_r1.fastq.gz"
    HOST_REMOVED_R2="${SAMPLE_ID}_host_removed_r2.fastq.gz"

    # Step 5: Align reads to the reference genome
    echo "Aligning reads for $SAMPLE_ID..."
    bwa mem -t $THREADS $REFERENCE $HOST_REMOVED_R1 $HOST_REMOVED_R2 > ${SAMPLE_ID}_aligned.sam

    # Step 6: Convert SAM to BAM, sort, and index
    echo "Converting SAM to sorted BAM for $SAMPLE_ID..."
    samtools view -Sb ${SAMPLE_ID}_aligned.sam > ${SAMPLE_ID}_aligned.bam
    samtools sort ${SAMPLE_ID}_aligned.bam -o ${SAMPLE_ID}_sorted.bam
    samtools index ${SAMPLE_ID}_sorted.bam

    # Clean up intermediate SAM file
    rm -f ${SAMPLE_ID}_aligned.sam

    # Step 7: Call variants
    echo "Calling variants for $SAMPLE_ID..."
    bcftools mpileup -Ou -f $REFERENCE ${SAMPLE_ID}_sorted.bam | bcftools call -mv -Oz -o ${SAMPLE_ID}_variants.vcf.gz
    bcftools index ${SAMPLE_ID}_variants.vcf.gz

    # Step 8: Generate consensus sequence
    echo "Generating consensus sequence for $SAMPLE_ID..."
    bcftools consensus -f $REFERENCE ${SAMPLE_ID}_variants.vcf.gz > ${SAMPLE_ID}_consensus.fasta

    # Step 9: Coverage analysis (Optional)
    echo "Calculating coverage for $SAMPLE_ID..."
    samtools depth ${SAMPLE_ID}_sorted.bam > ${SAMPLE_ID}_coverage.txt

    echo "Finished processing $SAMPLE_ID."
done < $SAMPLES_LIST

echo "All samples processed successfully."
