#!/bin/bash

module load R
module load samtools
export PATH=$HOME/bin/hmmcopy_utils/build/bin:$PATH

ID="4658_MM"
BAM="/scratch/dyerseq/hf180/4658_MM_work/${ID}_merged.bam"
OUTDIR="/scratch/dyerseq/hf180/4658_MM_work/results/downsampled_runs"
ICHOR_SCRIPT="/lustre/alice3/scratch/dyerseq/hf180/ichorCNA/scripts/runIchorCNA.R"
GCBASE="/scratch/dyerseq/hf180/ichorCNA/inst/extdata/gc_hg38_"
MAPBASE="/scratch/dyerseq/hf180/ichorCNA/inst/extdata/map_hg38_"
CENTROMERE="/scratch/dyerseq/hf180/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt"
NORMAL="c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)"
CHRS='c(1:22, "X")'
CHR_NORM_TRAIN='c(1:22)'
GENOME_BUILD="hg38"
GENOME_STYLE="UCSC"

mkdir -p "$OUTDIR"

for frac in 0.1 0.01 0.001; do
    frac_label=$(echo "$frac" | sed 's/0\.//')

    WORKDIR="$OUTDIR/ds${frac_label}"
    WIGDIR="$WORKDIR/wigs"
    RESULTDIR="$WORKDIR/results"
    mkdir -p "$WIGDIR" "$RESULTDIR"

    DSBAM="${WORKDIR}/${ID}_ds${frac_label}.bam"
    SORTED="${WORKDIR}/${ID}_ds${frac_label}.sorted.bam"

    samtools view -@ 8 -s "$frac" -b "$BAM" > "$DSBAM"
    samtools sort -@ 8 -o "$SORTED" "$DSBAM"
    samtools index "$SORTED"

    for binsize in 10000 50000 500000 1000000; do
        readCounter --window "$binsize" --quality 1 "$SORTED" > "$WIGDIR/${ID}.${binsize}.wig"
    done

    cd "$WIGDIR"
    for wig in *.wig; do
        binsize=$(echo "$wig" | sed -E 's/.*\.([0-9]+)\.wig/\1/')
        binsize_kb=$((binsize / 1000))
        outfile="${ID}.${binsize_kb}kb.filtered.wig"

        awk '
        /^fixedStep/ {
            split($0, a, "chrom=");
            split(a[2], b, " ");
            chrom = b[1];
            print_line = (chrom ~ /^chr([1-9][0-9]?|X)$/);
        }
        print_line { print }
        ' "$wig" > "$outfile"
    done

    for filtered in *.filtered.wig; do
        if [[ "$filtered" =~ \.([0-9]+)\.filtered\.wig$ ]]; then
            binsize="${BASH_REMATCH[1]}"
            binsize_kb=$((binsize / 1000))
            newname=$(echo "$filtered" | sed -E "s/\.${binsize}\.filtered\.wig$/.${binsize_kb}kb.filtered.wig/")
            mv "$filtered" "$newname"
        fi
    done

    for binsize in 10kb 50kb 500kb 1000kb; do
        WIG="${WIGDIR}/${ID}.${binsize}.filtered.wig"
        GCWIG="${GCBASE}${binsize}.wig"
        MAPWIG="${MAPBASE}${binsize}.wig"
        OUTBIN="${RESULTDIR}/bin${binsize}"
        mkdir -p "$OUTBIN"

        Rscript "$ICHOR_SCRIPT" \
            --id="${ID}_${frac_display}_downsampled_${binsize}" \
            --WIG="$WIG" \
            --gcWig="$GCWIG" \
            --mapWig="$MAPWIG" \
            --centromere="$CENTROMERE" \
            --normal="$NORMAL" \
            --genomeBuild="$GENOME_BUILD" \
            --genomeStyle="$GENOME_STYLE" \
            --outDir="$OUTBIN" \
            --chrs="$CHRS" \
            --chrNormalize="$CHR_NORM_TRAIN" \
            --chrTrain="$CHR_NORM_TRAIN"
    done
done

can you change it so that it does it for both MMS and RDV?? and also, change the files paths from MM to MMS
