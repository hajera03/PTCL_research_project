#!/bin/bash

#load required modules
module load samtools
export PATH=$PATH:/data/dyerseq/Nanopore_PTCL_data/minimap/minimap2

#define paths to viral reference genome
EBV="/lustre/alice3/scratch/dyerseq/hf180/sequences.fasta"
HTLV1="/lustre/alice3/scratch/dyerseq/hf180/HTLV1_sequence.fasta"

#indexing reference genomes if index files don't exist already
REFINDEX="/lustre/alice3/scratch/dyerseq/hf180/ref_indexed"
mkdir -p "$REFINDEX"
[ ! -f "$REFINDEX/EBV.mmi" ] && minimap2 -d "$REFINDEX/EBV.mmi" "$EBV"
[ ! -f "$REFINDEX/HTLV1.mmi" ] && minimap2 -d "$REFINDEX/HTLV1.mmi" "$HTLV1"

#sample name dictionary with corresponding fastq file locations
declare -A FASTQ_PATHS=(
  	["4408_RDV"]="/lustre/alice3/scratch/dyerseq/hf180/4408_RDV_work/fastq/combined.fastq"
  	["4658_MMS"]="/lustre/alice3/scratch/dyerseq/hf180/4658_MMS_work/fastq/combined.fastq"
)

#looping over each sample to perform alignment
for SAMPLE in "${!FASTQ_PATHS[@]}"; do
  	FASTQ="${FASTQ_PATHS[$SAMPLE]}"
  	OUTDIR="/lustre/alice3/scratch/dyerseq/hf180/${SAMPLE}_work/results/viral_detection"

  	mkdir -p "$OUTDIR"

  	#align to the human genome to remove host reads and then extract unmapped reads for viral alignment
   	#disabled to avoid missing weakly mapped viral reads
    
 	#minimap2 -ax sr map-ont "$REFINDEX/human.mmi" "$FASTQ" | samtools view -bS - > "$OUTDIR/${SAMPLE}_to_human.bam"
  		#samtools view -b -f 4 "$OUTDIR/${SAMPLE}_to_human.bam" > "$OUTDIR/${SAMPLE}_unmapped.bam"
  	#samtools fastq "$OUTDIR/${SAMPLE}_unmapped.bam" > "$OUTDIR/${SAMPLE}_unmapped.fastq"

  	#align to EBV
  	minimap2 -ax map-ont "$REFINDEX/EBV.mmi" "$FASTQ" | \
  		samtools view -bS - | samtools sort -o "$OUTDIR/${SAMPLE}_EBV.bam"
  	samtools index "$OUTDIR/${SAMPLE}_EBV.bam"

  	#align to HTLV1
	minimap2 -ax map-ont "$REFINDEX/HTLV1.mmi" "$FASTQ" | \
  		samtools view -bS - | samtools sort -o "$OUTDIR/${SAMPLE}_HTLV1.bam"
  	samtools index "$OUTDIR/${SAMPLE}_HTLV1.bam"
done

