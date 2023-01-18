#!/bin/bash
SRR=$1
THREADS=$2

prefetch -X 100000000 $SRR 
parallel-fastq-dump --split-files -t $THREADS -s $SRR --tmpdir /scratch1/bwu4/TEMP

bbduk.sh in1=${SRR}_1.fastq in2=${SRR}_2.fastq out1=${SRR}_1.clean.fastq out2=${SRR}_2.clean.fastq ref=TruSeq3-PE.fa ktrim=r k=23 \
        mink=8 hdist=1 qtrim=lr trimq=10 -minlen=80 tbo

bwa mem -t ${THREADS} CK2021.60.corrected.FOUR_TE_MASKED.fasta /zfs/socbd/bwu4/YU/VA_NGS/${SRR}_1.clean.fastq \
        /zfs/socbd/bwu4/YU/VA_NGS/${SRR}_2.clean.fastq | samtools sort -O BAM -o ${SRR}.TE_masked.srt.bam
