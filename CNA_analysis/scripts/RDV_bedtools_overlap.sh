#!/bin/bash

module load bedtools2

ID="4408_RDV"
PREFIX="RDV"
SEGDIR="/scratch/dyerseq/hf180/4408_RDV_work/results"
BED="/scratch/dyerseq/hf180/COSMIC/CGC.bed"
OUTDIR="/scratch/dyerseq/hf180/COSMIC/overlaps"

mkdir -p "$OUTDIR"

for binsize in 10kb 50kb 500kb 1000kb; do
    seg_file="${SEGDIR}/bin${binsize}/${ID}_${binsize}.seg"
    out_file="${OUTDIR}/${PREFIX}_${binsize}_overlaps.tsv"

    bedtools intersect \
        -a <(tail -n +2 "$seg_file" | awk -F'\t' '$5 != "NEUT"' | cut -f2-5) \
        -b "$BED" \
        -wa -wb \
        > "$out_file"
done
