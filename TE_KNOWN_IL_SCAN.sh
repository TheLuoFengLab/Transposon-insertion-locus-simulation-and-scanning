#!/bin/bash
SAMP=$1
IL=$2
NAME=$3
LEFT=$( awk -F '[:-]' '{print $2}' <<< "${IL}" )
RIGHT=$( awk -F '[:-]' '{print $3}' <<< "${IL}" )
samtools index ${SAMP}.TE_masked.srt.bam
samtools view -q 1 -h ${SAMP}.TE_masked.srt.bam ${IL} | awk -v x=${LEFT} -v y=${RIGHT} '
$17~/SA:Z:/ && ($4<x+150 || $4>y-150) { match($17,/SA:Z:([^,]*),([^,]*),([+-]),([^,]*),([^,]*),/,a);
if ($17~/,[0-9]+S[0-9]+M,/ && a[5]>0) {print a[1]"\t"a[2]-5"\t"a[2]+5"\t"$1};
if ($17~/,[0-9]+M[0-9]+S,/ && a[5]>0) {match($17,/,([0-9]+)M/,b); print a[1]"\t"a[2]+b[1]-5"\t"a[2]+b[1]+5"\t"$1};
if ($17~/,[0-9]+S[0-9]+M[0-9]+S,/ && a[5]>0) {match($17,/,([0-9]+)S([0-9]+)M([0-9]+)S,/,b); 
if (b[1]>b[3]) {print a[1]"\t"a[2]-5"\t"a[2]+5"\t"$1} else {print a[1]"\t"a[2]+b[2]-5"\t"a[2]+b[2]+5"\t"$1}}
}' > ${SAMP}_${NAME}.bed

bedtools sort -i ${SAMP}_${NAME}.bed | bedtools merge -d 100 -i - > ${SAMP}_${NAME}.merged.bed
