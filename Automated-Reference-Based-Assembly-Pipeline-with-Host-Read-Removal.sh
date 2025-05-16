#!/usr/bin/env bash
set -eo pipefail
trap 'echo "Error at line $LINENO"; exit 1' ERR

# Configuration
THREADS=$(nproc)
REFERENCE="reference.fasta"
HOST_REFERENCE="host_reference.fasta"
SAMPLES_LIST="samples.txt"
OUTDIR="analysis_results"
QUALITY=20       # Minimum quality score for trimming
MIN_LENGTH=50    # Minimum read length after trimming

# Initialize directories
mkdir -p "${OUTDIR}/"{qc,trimmed,host_filtered,alignment,variants,consensus,coverage}

# Dependency check
declare -A TOOLS=(
    [bwa]="conda install -c bioconda bwa"
    [minimap2]="conda install -c bioconda minimap2"
    [fastqc]="conda install -c bioconda fastqc"
    [fastp]="conda install -c bioconda fastp"
    [samtools]="conda install -c bioconda samtools"
    [bcftools]="conda install -c bioconda bcftools"
    [hostile]="pip install hostile"
)

for tool in "${!TOOLS[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "ERROR: $tool not found. Install with: ${TOOLS[$tool]}"
        exit 1
    fi
done

# Index references
index_references() {
    echo "Indexing references..."
    [ -f "${REFERENCE}.bwt" ] || bwa index "$REFERENCE"
    [ -f "reference.mmi" ] || minimap2 -d reference.mmi "$REFERENCE"
    [ -f "${HOST_REFERENCE}.bwt" ] || bwa index "$HOST_REFERENCE"
}

process_sample() {
    local sample=$1
    echo "Processing ${sample}..."
    
    # Input files
    local raw_r1="${sample}_R1.fastq.gz"
    local raw_r2="${sample}_R2.fastq.gz"
    
    # Quality control
    fastqc -o "${OUTDIR}/qc" "$raw_r1" "$raw_r2"
    
    # Adapter trimming
    fastp -w "$THREADS" \
        -i "$raw_r1" -I "$raw_r2" \
        -o "${OUTDIR}/trimmed/${sample}_trimmed_R1.fastq.gz" \
        -O "${OUTDIR}/trimmed/${sample}_trimmed_R2.fastq.gz" \
        -q $QUALITY -l $MIN_LENGTH \
        --json "${OUTDIR}/trimmed/${sample}_fastp.json" \
        --html "${OUTDIR}/trimmed/${sample}_fastp.html"

    # Host read removal
    hostile --host-index "$HOST_REFERENCE" \
        --threads "$THREADS" \
        --r1 "${OUTDIR}/trimmed/${sample}_trimmed_R1.fastq.gz" \
        --r2 "${OUTDIR}/trimmed/${sample}_trimmed_R2.fastq.gz" \
        --output-prefix "${OUTDIR}/host_filtered/${sample}_host_removed"

    # Alignment and sorting
    bwa mem -t "$THREADS" "$REFERENCE" \
        "${OUTDIR}/host_filtered/${sample}_host_removed_r1.fastq.gz" \
        "${OUTDIR}/host_filtered/${sample}_host_removed_r2.fastq.gz" |
    samtools sort -@ "$THREADS" -o "${OUTDIR}/alignment/${sample}_sorted.bam"
    
    samtools index "${OUTDIR}/alignment/${sample}_sorted.bam"

    # Variant calling
    bcftools mpileup -Ou -f "$REFERENCE" \
        "${OUTDIR}/alignment/${sample}_sorted.bam" |
    bcftools call -mv -Oz \
        -o "${OUTDIR}/variants/${sample}_variants.vcf.gz"
    
    bcftools index "${OUTDIR}/variants/${sample}_variants.vcf.gz"

    # Consensus generation
    bcftools consensus -f "$REFERENCE" \
        "${OUTDIR}/variants/${sample}_variants.vcf.gz" \
        > "${OUTDIR}/consensus/${sample}_consensus.fasta"

    # Coverage analysis
    samtools depth -aa "${OUTDIR}/alignment/${sample}_sorted.bam" |
    awk '{sum+=$3} END {print "Average coverage:", sum/NR}' \
        > "${OUTDIR}/coverage/${sample}_coverage.txt"
}

# Main execution
index_references
while read -r sample; do
    process_sample "$sample"
done < "$SAMPLES_LIST"

echo "Pipeline completed successfully. Results in: ${OUTDIR}"
