#!/bin/bash
#SBATCH --job-name=merge_index_bams         # Job name for SLURM queue
#SBATCH --nodes=1                           # Run on a single node
#SBATCH --ntasks=1                          # Number of tasks
#SBATCH --cpus-per-task=16                  # Use 16 threads
#SBATCH --mem=200G                          # Memory allocation
#SBATCH --time=10:00:00                     # Max walltime (hh:mm:ss)
#SBATCH --mail-type=BEGIN,END,FAIL          # Get emails on job begin, end, and fail
#SBATCH --mail-user=hf180@student.le.ac.uk  # Your email

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
