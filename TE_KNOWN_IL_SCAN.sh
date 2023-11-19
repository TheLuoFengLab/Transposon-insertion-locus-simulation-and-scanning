#!/bin/bash
# Shebang line to indicate the script should be run using Bash

SAMP=$1
# Assigns the first command-line argument to the variable SAMP, typically a sample identifier

IL=$2
# Assigns the second command-line argument to the variable IL, representing an interval in a specific format

NAME=$3
# Assigns the third command-line argument to the variable NAME, used as a tag for naming output files

LEFT=$( awk -F '[:-]' '{print $2}' <<< "${IL}" )
RIGHT=$( awk -F '[:-]' '{print $3}' <<< "${IL}" )
# Extracts the left and right boundaries from the IL variable using awk

samtools index ${SAMP}.TE_masked.srt.bam
# Indexes the sorted BAM file using samtools for efficient access

samtools view -q 1 -h ${SAMP}.TE_masked.srt.bam ${IL} | \
awk -v x=${LEFT} -v y=${RIGHT} '$17~/SA:Z:/ && ($4<x+150 || $4>y-150) { match($17,/SA:Z:([^,]*),([^,]*),([+-]),([^,]*),([^,]*),/,a); 
if ($17~/,[0-9]+S[0-9]+M,/ && a[5]>0) {print a[1]"\t"a[2]-5"\t"a[2]+5"\t"$1}; 
if ($17~/,[0-9]+M[0-9]+S,/ && a[5]>0) {match($17,/,([0-9]+)M/,b); print a[1]"\t"a[2]+b[1]-5"\t"a[2]+b[1]+5"\t"$1}; 
if ($17~/,[0-9]+S[0-9]+M[0-9]+S,/ && a[5]>0) {match($17,/,([0-9]+)S([0-9]+)M([0-9]+)S,/,b);  
if (b[1]>b[3]) {print a[1]"\t"a[2]-5"\t"a[2]+5"\t"$1} else {print a[1]"\t"a[2]+b[2]-5"\t"a[2]+b[2]+5"\t"$1}}}' \
> ${SAMP}_${NAME}.bed
# Line 1: Uses 'samtools view' to extract reads from the BAM file within the specified interval (IL)
# with a minimum mapping quality of 1 (-q 1). The -h flag includes the header in the output.
# The output is piped into an awk command for further processing.
# Line 2: Begins an awk command, passing in variables x and y set to the values of LEFT and RIGHT.
# Filters reads that contain supplementary alignments (indicated by "SA:Z:" in the 17th field)
# and whose start position ($4) is within 150 bases of the interval edges (LEFT and RIGHT).
# Uses the match function to extract supplementary alignment information from the 17th field.
# Captured groups are stored in the array 'a'.
# Line 3: Checks if the CIGAR string in the supplementary alignment indicates soft clipping (S) followed by matching (M),
# and if the mapping quality (stored in a[5]) is greater than 0. Prints the chromosome (a[1]), start position minus 5 (a[2]-5), end position plus 5 (a[2]+5),
# and the read name ($1), separated by tabs.
# Line 4: Matches a pattern where matching (M) is followed by soft clipping (S) in the CIGAR string.
# Captures the length of the matching segment.
# Calculates and prints the appropriate genomic coordinates and read name.
# Line 5: Matches a CIGAR string pattern with soft clipping, matching, and then soft clipping.
# Captures the lengths of these segments.
# Line 6: Determines which clipping is larger and prints coordinates and read name accordingly.
# Line 7: The output from awk is redirected to a BED file named based on the sample identifier (SAMP) and the provided name tag (NAME).

bedtools sort -i ${SAMP}_${NAME}.bed | bedtools merge -d 100 -i - > ${SAMP}_${NAME}.merged.bed
# This command first sorts the BED file using bedtools sort
# Then, it merges intervals within 100 bases of each other using bedtools merge
# The final merged BED file is saved with the sample identifier and name tag
