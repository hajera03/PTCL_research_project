#!/bin/bash

module load R
export PATH=$HOME/bin/hmmcopy_utils/build/bin:$PATH

ID="4408_RDV"
BAM="/scratch/dyerseq/hf180/4408_RDV_work/${ID}_merged.bam"
WIGDIR="/scratch/dyerseq/hf180/4408_RDV_work/wigs"
RESULTDIR="/scratch/dyerseq/hf180/4408_RDV_work/results"
ICHOR_SCRIPT="/lustre/alice3/scratch/dyerseq/hf180/ichorCNA/scripts/runIchorCNA.R"
GCBASE="/scratch/dyerseq/hf180/ichorCNA/inst/extdata/gc_hg38_"
MAPBASE="/scratch/dyerseq/hf180/ichorCNA/inst/extdata/map_hg38_"
CENTROMERE="/scratch/dyerseq/hf180/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt"
#COVERAGE=13.2112
NORMAL="c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)"
#SCSTATES="c(1,3)"
CHRS='c(1:22, "X")'
CHR_NORM_TRAIN='c(1:22)'
GENOME_BUILD="hg38"
GENOME_STYLE="UCSC"

for binsize in 10000 50000 500000 1000000; do
    WIG="${WIGDIR}/${ID}.${binsize}.wig"
    if [[ ! -f "$WIG" ]]; then
        readCounter --window "$binsize" --quality 1 "$BAM" > "$WIG"
    else
        echo "Skipping readCounter: $WIG already exists"
    fi
done

cd "$WIGDIR"
for wig in *.wig; do
    binsize=$(echo "$wig" | sed -E 's/.*\.([0-9]+)\.wig/\1/')
    binsize_kb=$((binsize / 1000))kb
    outfile="${ID}.${binsize}.filtered.wig"

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
        binsize_kb="$((binsize / 1000))kb"
        newname=$(echo "$filtered" | sed -E "s/\.${binsize}\.filtered\.wig$/.${binsize_kb}.filtered.wig/")
        mv "$filtered" "$newname"
    fi
done

# Step 4: Run ichorCNA for each bin size
for binsize in 10kb 50kb 500kb 1000kb; do
    WIG="${WIGDIR}/${ID}.${binsize}.filtered.wig"
    GCWIG="${GCBASE}${binsize}.wig"
    MAPWIG="${MAPBASE}${binsize}.wig"
    OUTDIR="${RESULTDIR}/bin${binsize}"
    mkdir -p "$OUTDIR"

    Rscript "$ICHOR_SCRIPT" \
        --id="${ID}_${binsize}" \
        --WIG="$WIG" \
        --gcWig="$GCWIG" \
        --mapWig="$MAPWIG" \
        --centromere="$CENTROMERE" \
        --normal="$NORMAL" \
        --genomeBuild="$GENOME_BUILD" \
        --genomeStyle="$GENOME_STYLE" \
        --outDir="$OUTDIR" \
        --chrs="$CHRS" \
        --chrNormalize="$CHR_NORM_TRAIN" \
        --chrTrain="$CHR_NORM_TRAIN"
done
