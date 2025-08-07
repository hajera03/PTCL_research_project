#!/bin/bash
module load kraken2

#define paths
DB_PATH="/lustre/alice3/scratch/dyerseq/hf180/minikraken_8GB_20200312"
OUTDIR="/lustre/alice3/scratch/dyerseq/hf180/kraken2_results"
mkdir -p "$OUTDIR"
samples=("4408_RDV" "4658_MMS")

#loop through each sample
for SAMPLE in "${samples[@]}"; do
    FASTQ="/lustre/alice3/scratch/dyerseq/hf180/${SAMPLE}_work/fastq/combined.fastq.gz"
    REPORT="${OUTDIR}/${SAMPLE}_full_report.txt"
    OUTPUT="${OUTDIR}/${SAMPLE}_full_output.txt"
    CLASSIFIED="${OUTDIR}/${SAMPLE}_full_classified.fastq"

    kraken2 \
        --db "$DB_PATH" \
        --threads 16 \
        --report "$REPORT" \
        --output "$OUTPUT" \
        --classified-out "$CLASSIFIED" \
        "$FASTQ"
done
