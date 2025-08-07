#this script merges and indexes multiple bam files for two ctDNA samples : 4408_RDV and 4658_MMS using samtools
#!/bin/bash

#load necessary modules
module load samtools

#defining sample IDs
SAMPLES=("4408_RDV" "4658_MMS")

#loop over each sample
for SAMPLE in "${SAMPLES[@]}"; do
    #defining working directory
    workdir="/scratch/dyerseq/hf180/${SAMPLE}_work"
    cd "$workdir"

    #create a list of all bam files
    ls *.bam > bamlist.txt

    #output filename
    merged="${SAMPLE}_merged.bam"

    #merge all bam files in bamlist.txt using 16 threads
    samtools merge -@ 16 -b bamlist.txt "$merged"

    #check the validity of the merged bam file before indexing
    if samtools quickcheck "$merged"; then
        #creating bam index (.bai)
        samtools index -@ 16 "$merged"
    fi
done
