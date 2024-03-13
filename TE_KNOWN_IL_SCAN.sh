#!/bin/bash
# Shebang line to indicate the script should be run using Bash

SAMP=$1
# Assigns the first command-line argument to the variable SAMP, typically a sample identifier

IL=$2
# IL format should be chromosome:start-end, for instance chr2:10000-10009, TE insertion usually renders 
# duplication of a short sequence (<=20 bp) at the insertion locus, therefore end - start should be less than 20

NAME=$3
# Assigns the third command-line argument to the variable NAME, used as a tag for naming output files

READL=$4
# Read length

TE=$5
# ID of the TE in the masked reference fasta file 

LEFT=$( awk -F '[:-]' '{print $2}' <<< "${IL}" )
RIGHT=$( awk -F '[:-]' '{print $3}' <<< "${IL}" )
# Extracts the left and right boundaries from the IL variable using awk

samtools index ${SAMP}.TE_masked.srt.bam
# Indexes the sorted BAM file using samtools for efficient access

# Ourput records supporting 
samtools view -q 1 -h ${SAMP}.TE_masked.srt.bam ${IL} | \
awk -v x=$LEFT -v y=$RIGHT -v z=$READL '$0~/SA:Z:/ && ($4<x+z || $4>y-z) { 
    match($0,/SA:Z:([^,]*),([^,]*),([+-]),([^,]*),([^,]*),/,a); 
    if (a[4]~/^[0-9]+S[0-9]+M$/ && a[5]>0) {
        print a[1]"\t"a[2]-5"\t"a[2]+5"\t"$1
    } else if (a[4]~/^[0-9]+M[0-9]+S$/ && a[5]>0) {
        match(a[4],/^([0-9]+)M/,b); 
        print a[1]"\t"a[2]+b[1]-5"\t"a[2]+b[1]+5"\t"$1
    } else if (a[4]~/^[0-9]+S[0-9]+M[0-9]+S$/ && a[5]>0) {
        match(a[4],/^([0-9]+)S([0-9]+)M([0-9]+)S$/,b);  
        if (b[1]>b[3]) {
            print a[1]"\t"a[2]-5"\t"a[2]+5"\t"$1
        } else {
            print a[1]"\t"a[2]+b[2]-5"\t"a[2]+b[2]+5"\t"$1
        }
    }
}' > ${SAMP}_${NAME}.bed


