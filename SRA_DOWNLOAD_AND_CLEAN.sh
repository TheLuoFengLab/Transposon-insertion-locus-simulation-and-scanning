#!/bin/bash
# Shebang line to indicate the script should be run using Bash

SRR=$1
# Assign the first command-line argument to the variable SRR (usually a sample identifier)

THREADS=$2
# Assign the second command-line argument to the variable THREADS (number of CPU threads to use)

prefetch -X 100000000 $SRR 
# Use the 'prefetch' command from the SRA Toolkit to download data for the SRR ID with a cache size of 100,000,000 bytes

parallel-fastq-dump --split-files -t $THREADS -s $SRR --tmpdir /scratch1/bwu4/TEMP
# Use 'parallel-fastq-dump' to convert SRA files to FASTQ format, splitting them into separate files (if paired-end), 
# using specified number of threads, and using a temporary directory for intermediate files

bbduk.sh in1=${SRR}_1.fastq in2=${SRR}_2.fastq out1=${SRR}_1.clean.fastq out2=${SRR}_2.clean.fastq ref=TruSeq3-PE.fa ktrim=r k=23 \
        mink=8 hdist=1 qtrim=lr trimq=10 -minlen=80 tbo
# Run 'bbduk.sh' (from BBTools) for quality control and trimming. It takes input FASTQ files (in1 and in2), 
# outputs cleaned FASTQ files (out1 and out2), trims adapters based on the provided reference (TruSeq3-PE.fa),
# applies k-mer trimming, sets quality trimming parameters, and filters out reads below 80bp after trimming

bwa mem -t ${THREADS} CK2021.60.corrected.FOUR_TE_MASKED.fasta /zfs/socbd/bwu4/YU/VA_NGS/${SRR}_1.clean.fastq \
        /zfs/socbd/bwu4/YU/VA_NGS/${SRR}_2.clean.fastq | samtools sort -O BAM -o ${SRR}.TE_masked.srt.bam
# Execute 'bwa mem' for sequence alignment using the specified number of threads, with the provided reference genome (CK2021.60.corrected.FOUR_TE_MASKED.fasta),
# and the cleaned FASTQ files as input. The output is piped to 'samtools sort', which sorts alignments and outputs a sorted BAM file named ${SRR}.TE_masked.srt.bam

