# Transposon insertion locus (IL) scanning
## Software installation



# # File content for FOUR_TE_INSERTION.sh (github)

```bash
#!/bin/bash
SAMP=$1
TE=$2
NAME=$3
LEFT=$( awk -F '[:-]' '{print $2}' <<< "${TE}" )
RIGHT=$( awk -F '[:-]' '{print $3}' <<< "${TE}" )
samtools index ${SAMP}.TE_masked.srt.bam
samtools view -q 1 -h ${SAMP}.TE_masked.srt.bam ${TE} | awk -v x=${LEFT} -v y=${RIGHT} '
$17~/SA:Z:/ && ($4<x+150 || $4>y-150) { match($17,/SA:Z:([^,]*),([^,]*),([+-]),([^,]*),([^,]*),/,a);
if ($17~/,[0-9]+S[0-9]+M,/ && a[5]>0) {print a[1]"\t"a[2]-5"\t"a[2]+5"\t"$1};
if ($17~/,[0-9]+M[0-9]+S,/ && a[5]>0) {match($17,/,([0-9]+)M/,b); print a[1]"\t"a[2]+b[1]-5"\t"a[2]+b[1]+5"\t"$1};
if ($17~/,[0-9]+S[0-9]+M[0-9]+S,/ && a[5]>0) {match($17,/,([0-9]+)S([0-9]+)M([0-9]+)S,/,b); 
if (b[1]>b[3]) {print a[1]"\t"a[2]-5"\t"a[2]+5"\t"$1} else {print a[1]"\t"a[2]+b[2]-5"\t"a[2]+b[2]+5"\t"$1}}
}' > ${SAMP}_${NAME}.bed

bedtools sort -i ${SAMP}_${NAME}.bed | bedtools merge -d 100 -i - > ${SAMP}_${NAME}.merged.bed
```

# # Generating PBS scripts for the rest sweet orange SRRs

```bash
for SRR in $( cat REST.SRR | grep -v 1025 ); do
echo "#PBS -N $SRR
#PBS -l select=1:ncpus=16:mem=60gb:interconnect=1g,walltime=72:00:00
#PBS -j oe
source ~/.bashrc
conda activate ecDNA
cd /zfs/socbd/bwu4/YU/annotation/FOUR_TE
./PRE_CLE.sh ${SRR} 16" > ${SRR}.pbs
done
```

# # Scanning TE insertion loci in the rest SRR bam

```bash
for SAM in $( cat REST.SRR ); do
./FOUR_TE_INSERTION.sh $SAM chr1A:16314816-16321163 MU1
done

for SAM in $( cat REST.SRR ); do
./FOUR_TE_INSERTION.sh $SAM chr5B:2319829-2326754 MU2
done

for SAM in $( cat FASTQ.list ); do
for TE in MU1 MU2 ; do
if test ! -f ${SAM}_${TE}.merged.filt.bed ; then
    bedtools intersect -wa -a ${SAM}_${TE}.merged.bed -b ${SAM}_${TE}.bed | uniq -c | \
        sed "s/ \+/\t/g;s/^\t//g" | awk '$1>1' | cut -f 2,3,4 > ${SAM}_${TE}.merged.filt.bed
fi
done
done

cat *_MU1.merged.filt.bed | bedtools sort -i - | bedtools merge -d 100 -i - > ALL_MU1.merged.filt.merged.bed
cat *_MU2.merged.filt.bed | bedtools sort -i - | bedtools merge -d 100 -i - > ALL_MU2.merged.filt.merged.bed

for SAM in $( cat FASTQ.list ); do
bedtools intersect -c -a ALL_MU1.merged.filt.merged.bed -b ${SAM}_MU1.merged.filt.bed > ${SAM}_MU1.count.tsv
bedtools intersect -c -a ALL_MU2.merged.filt.merged.bed -b ${SAM}_MU2.merged.filt.bed > ${SAM}_MU2.count.tsv
done

cut -f 1,2,3 CK_MU1.count.tsv > ALL_MU1.count.tsv
cut -f 1,2,3 CK_MU2.count.tsv > ALL_MU2.count.tsv

for SAM in $( cat FASTQ.list ); do
cut -f 4 ${SAM}_MU1.count.tsv | paste ALL_MU1.count.tsv - > ALL_MU1.count.tsv.temp
mv ALL_MU1.count.tsv.temp ALL_MU1.count.tsv
cut -f 4 ${SAM}_MU2.count.tsv | paste ALL_MU2.count.tsv - > ALL_MU2.count.tsv.temp
mv ALL_MU2.count.tsv.temp ALL_MU2.count.tsv
done

sed ':a;N;$!ba;s/\n/\t/g' FASTQ.list > header.txt  # replace \n with \t
vi header.txt # add three columns to the header "CHR    START   END"
cat header.txt ALL_MU1.count.tsv > ALL_MU1.count.tsv.temp
mv ALL_MU1.count.tsv.temp ALL_MU1.count.tsv
cat header.txt ALL_MU2.count.tsv > ALL_MU2.count.tsv.temp
mv ALL_MU2.count.tsv.temp ALL_MU2.count.tsv

awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' CK_MU1.count.tsv > MU1.INS_LOCI.bed
awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' CK_MU2.count.tsv > MU2.INS_LOCI.bed

bedtools intersect -wao -a MU2.INS_LOCI.bed -b ALL_FOUR_TE.EXT100.bed | awk '$5~/chr/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5":"$6"-"$7} $5!~/chr/ {print $1"\t"$2"\t"$3"\t"$4"\t"$1":"$2"-"$3}' > MU2.ID_CONVERSION.tsv
bedtools intersect -wao -a MU1.INS_LOCI.bed -b ALL_FOUR_TE.EXT100.bed | awk '$5~/chr/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5":"$6"-"$7} $5!~/chr/ {print $1"\t"$2"\t"$3"\t"$4"\t"$1":"$2"-"$3}' > MU1.ID_CONVERSION.tsv

cut -f 5 MU1.ID_CONVERSION.tsv | sort -u > MU1_INS_LOCI_FINAL.tsv
cut -f 5 MU2.ID_CONVERSION.tsv | sort -u > MU2_INS_LOCI_FINAL.tsv

awk '{match($1,/([^:]*):([^-]*)-([^-]*)/,a);print a[1]"\t"a[2]"\t"a[3]"\t"$1}' MU1_INS_LOCI_FINAL.tsv > MU1_INS_LOCI_FINAL.bed
awk '{match($1,/([^:]*):([^-]*)-([^-]*)/,a);print a[1]"\t"a[2]"\t"a[3]"\t"$1}' MU2_INS_LOCI_FINAL.tsv > MU2_INS_LOCI_FINAL.bed

bedtools intersect -wo -a MU1_INS_LOCI_FINAL.bed -b CK2021.FINAL.0403.gene.bed | awk '{print $4"\t"$8}' | sort -u > MU1.INS_LOCI.GENE.tsv
bedtools intersect -wo -a MU2_INS_LOCI_FINAL.bed -b CK2021.FINAL.0403.gene.bed | awk '{print $4"\t"$8}' | sort -u > MU2.INS_LOCI.GENE.tsv

bedtools intersect -wo -a MU1_INS_LOCI_FINAL.bed -b CK2021.FINAL.0403.CDS.bed | awk '{print $4"\t"$8}' | sort -u > MU1.INS_LOCI.CDS.tsv
bedtools intersect -wo -a MU2_INS_LOCI_FINAL.bed -b CK2021.FINAL.0403.CDS.bed | awk '{print $4"\t"$8}' | sort -u > MU2.INS_LOCI.CDS.tsv

bedtools intersect -wo -a MU1_INS_LOCI_FINAL.bed -b CK2021.FINAL.0403.UTR.bed  | awk '{print $4"\t"$8}' | sort -u > MU1.INS_LOCI.UTR.tsv
bedtools intersect -wo -a MU2_INS_LOCI_FINAL.bed -b CK2021.FINAL.0403.UTR.bed  | awk '{print $4"\t"$8}' | sort -u > MU2.INS_LOCI.UTR.tsv

bedtools intersect -wo -a MU1_INS_LOCI_FINAL.bed -b CK2021.FINAL.0403.gene.upstream.bed | awk '{print $4"\t"$8}' | sort -u > MU1.INS_LOCI.UPSTREAM.tsv
bedtools intersect -wo -a MU2_INS_LOCI_FINAL.bed -b CK2021.FINAL.0403.gene.upstream.bed | awk '{print $4"\t"$8}' | sort -u > MU2.INS_LOCI.UPSTREAM.tsv
```

# # Sequencing depth analysis

```bash
for SRR in $( cat FASTQ.list ); do
samtools view -@ $THREADS -b -q 1 -o ${SRR}.TE_masked.srt.q1.bam ${SRR}.TE_masked.srt.q1.bam
bedtools coverage -mean -a CK2021.100_2K.srt.wins -b ${SRR}.TE_masked.srt.q1.bam > ${SRR}.100_2K.cov
done

sed ':a;N;$!ba;s/\n/\t/g' FASTQ.list | echo "CHR\tSTART\tEND\tNAME\t"-
	
```

# # Identification of Mu1-3 in different citrus assemblies

```bash
ls /zfs/socbd/bwu4/REF/
# CCL  CSI  HKC_Box  HWB  KUMQUAT  MSYJ  PTR  XJC_Ci  XZ_Citron  ZK_Pt

cd /zfs/socbd/bwu4/YU/annotation/FOUR_TE/TE_IN_CITRUS_GENOMES
makeblastdb -in /zfs/socbd/bwu4/REF/KUMQUAT/GJ.contig.fa -dbtype nucl -parse_seqids -out GJ
mv GJ* /zfs/socbd/bwu4/REF/KUMQUAT/
blastn -query ../FOUR_TE.fasta -db /zfs/socbd/bwu4/REF/KUMQUAT/GJ -num_threads 40 -evalue 1E-3 -perc_identity 0.9 -outfmt 6 > FOUR_TE_GJ.results

makeblastdb -in /zfs/socbd/bwu4/REF/HKC_Box/HKC.scaffold.fa -dbtype nucl -parse_seqids -out HKC
mv HKC* /zfs/socbd/bwu4/REF/HKC_Box/
blastn -query ../FOUR_TE.fasta -db /zfs/socbd/bwu4/REF/HKC_Box/HKC -num_threads 40 -evalue 1E-3 -perc_identity 0.9 -outfmt 6 > FOUR_TE_HKC.results

makeblastdb -in /zfs/socbd/bwu4/REF/XJC_Ci/XJC.scaffold.fa -dbtype nucl -parse_seqids -out XJC
mv XJC* /zfs/socbd/bwu4/REF/XJC_Ci/
blastn -query ../FOUR_TE.fasta -db /zfs/socbd/bwu4/REF/XJC_Ci/XJC -num_threads 40 -evalue 1E-3 -perc_identity 0.9 -outfmt 6 > FOUR_TE_XJC.results

makeblastdb -in /zfs/socbd/bwu4/REF/XZ_Citron/XZ.scaffold.fa -dbtype nucl -parse_seqids -out XZ
mv XZ* /zfs/socbd/bwu4/REF/XZ_Citron/
blastn -query ../FOUR_TE.fasta -db /zfs/socbd/bwu4/REF/XZ_Citron/XZ -num_threads 40 -evalue 1E-3 -perc_identity 0.9 -outfmt 6 > FOUR_TE_XZ.results

for RESU in FOUR_TE_GJ.results FOUR_TE_HKC.results FOUR_TE_XJC.results FOUR_TE_XZ.results ; do
awk '$1=="LINE1" && ($7<50 || $8>6300) {print}' $RESU > ${RESU}.Mu1_TEMP
awk '$1=="LINE2" && ($7<50 || $8>6900) {print}' $RESU > ${RESU}.Mu2_TEMP
done

bedtools getfasta -s -nameOnly -fi /zfs/socbd/bwu4/REF/XJC_Ci/XJC.scaffold.fa -fo XJC_4TE.fasta -bed XJC_4TE.bed  

bedtools getfasta -s -nameOnly -fi /zfs/socbd/bwu4/REF/XZ_Citron/XZ.scaffold.fa -fo XZ_4TE.fasta -bed XZ_4TE.bed

bedtools getfasta -s -nameOnly -fi /zfs/socbd/bwu4/REF/KUMQUAT/GJ.contig.fa -fo HK_4TE.fasta -bed HK_4TE.bed

bedtools getfasta -s -nameOnly -fi /zfs/socbd/bwu4/REF/HKC_Box/HKC.scaffold.fa -fo HKC_4TE.fasta -bed HKC_4TE.bed
```

# # Mask reference genome

```bash
mkdir /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW 
cd /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW 
infoseq -only -name -length CK2021.60.corrected.A.fasta > CK2021.A.txt
infoseq -only -name -length CK2021.60.corrected.B.fasta > CK2021.B.txt
sed -i 1d CK2021.A.txt
sed -i 1d CK2021.B.txt
sed -i 's/  / /g' *.txt # Replace double spaces with single space multiple times
sed -i "s/ /\t/g" *.txt
bedtools makewindows -g CK2021.A.txt -w 50 -s 10 > CK2021.A.50_10.wins
bedtools getfasta -fi CK2021.60.corrected.A.fasta -bed CK2021.A.50_10.wins -fo CK2021.A.50_10.fasta

bwa mem -t 8 CK2021.60.corrected.B.fasta CK2021.A.50_10.fasta > CK2021.A.50_10.B.sam
awk '{if ($6~/50M/ && $12~/NM:i:0/) print $3"\t"$4-1"\t"$4+49}'  CK2021.A.50_10.B.sam > CK2021.B.DUP.bed
bedtools sort -i CK2021.B.DUP.bed | bedtools merge -i - > CK2021.B.DUP.merged.bed
awk '$3-$2>=101 {print $1"\t"$2+50"\t"$3-50}' CK2021.B.DUP.merged.bed > CK2021.B.TO_MASK.bed
cut -f 1,2,3 CORRECTED_DVS_4TE.bed > CORRECTED_DVS_4TE.1-3.bed
# add the 5 non-intact MULE regions in CORRECTED_DVS_4TE.1-3.bed
cat CORRECTED_DVS_4TE.1-3.bed CK2021.B.TO_MASK.bed | bedtools sort -i - | bedtools merge -i - > CK2021.ALLELE_Mu.MASKED.bed
bedtools maskfasta -fi CK2021.ALLELE_Mu.MASKED.fasta -bed CK2021.B.TO_MASK.bed -fo CK2021.09282021.MASKED.fasta

# The final masked regions are recorded in the file CK2021.ALLELE_MU.MASKED.bed
```

# # Mapping cleaned fastq to masked reference with separate mutators

```bash
ln /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/CK2021.ALLELE_Mu.MASKED.fasta
bwa index CK2021.ALLELE_Mu.MASKED.fasta
bwa mem -t 18 CK2021.09282021.MASKED.fasta /zfs/socbd/bwu4/YU/VA_NGS/${SAMP}_1.clean.fastq \
        /zfs/socbd/bwu4/YU/VA_NGS/${SAMP}_2.clean.fastq | samtools sort -@ 2 -O BAM -o ${SAMP}.TE_masked.srt.bam
```

# # Filter_Split_Reads.py

```python
#!/usr/bin/env python3

import pandas as pd
import sys
# sys.argv[1] : list of informative split mapped reads
# sys.argv[2] : path of two-column SAM file

SPLIT = pd.read_table(sys.argv[1], names=['ReadID'])
SAM = pd.read_table(sys.argv[2],names=['ReadID','Merged_columns'])
SAM.loc[~SAM['ReadID'].isin(SPLIT['ReadID'])].to_csv('SPLIT_FILT.'+sys.argv[2],sep="\t", index=False, header=False)
```

# # Content of TE_INSERT_READ.sh

```bash
#!/bin/bash

# Read transposable element (TE) end information from file $1 into the array a, with each line as an element
# The file format of the first argument (**space as seperator**; n denotes the full length of the TE in the reference): 
   # TE_ID 1 LEFT
   # TE_ID n RIGHT   
readarray -t a < $1

# SRR ID or sample ID, the sorted bam should be named as ${SAMP}.TE_masked.srt.bam
SAMP=$2
# Number of threads, must be integer
TR=$3
# Single-end sequencing read length 
READ_LEN=$4
# Average sequencing library insertion segment length
INS_LEN=$5
MAX_INS=$((INS_LEN*2))

# Common text in TE names. For instance, 'Mu' in TEs with names 'Mu1_1, Mu1_2, Mu1_3, Mu2_1, Mu2_2' in the reference fasta file
TE=$6

# Reference fasta file including modified Mu sequences
REF=$7

# Initialization of variable i
i=0

# Mapping cleaned fastq to the reference if the sorted bam does not exist
if [[ ! -f "${SAMP}.TE_masked.srt.bam" ]]; then
echo "Start mapping ${SAMP}_1.clean.fastq and ${SAMP}_1.clean.fastq to the masked reference"
	# Index reference fasta file if the index file does not exist
	if [[ ! -f ${REF} ]]; then
	echo "Error : neither the bam file of ${SAMP} nor the reference fasta file exist."	
	exit 1
	fi
	if [[ ! -f "${REF}.sa" ]]; then
		bwa index $REF
	fi

	if [[ -f "${SAMP}_1.clean.fastq" ]] && [[ -f "${SAMP}_2.clean.fastq" ]] ; then
		bwa mem -t $TR ${REF} ${SAMP}_1.clean.fastq \
		        ${SAMP}_2.clean.fastq | samtools sort -@ $TR -O BAM -o ${SAMP}.TE_masked.srt.bam
	elif [[ -f "${SAMP}_1.clean.fastq.gz" ]] && [[ -f "${SAMP}_2.clean.fastq.gz" ]] ; then
		bwa mem -t $TR ${REF} ${SAMP}_1.clean.fastq.gz \
		        ${SAMP}_2.clean.fastq.gz | samtools sort -@ $TR -O BAM -o ${SAMP}.TE_masked.srt.bam
	else
	echo "Error : no bam file or fastq file of ${SAMP} exist in current working directory."
	echo "Either ${SAMP}_[12].clean.fastq.gz or ${SAMP}_[12].clean.fastq files should be supplied"
	exit 2
	fi
fi

#Index the sorted BAM file
if [[ ! -f "${SAMP}.TE_masked.srt.bam.bai" ]]; then
samtools index -@ $TR ${SAMP}.TE_masked.srt.bam
fi

# Delete previously produced files
rm -f ${SAMP}.SPLIT_SUP.bed
rm -f ${SAMP}.DISC_SUP.bed
rm -f ${SAMP}.*.sam

# Re-sort the bam via read name
samtools view -h -f 1 -F 2318 ${SAMP}.TE_masked.srt.bam | samtools sort -@ $((TR-1)) -n -o ${SAMP}.disc.readname_srt.bam -

while [[ $i -lt ${#a[@]} ]]; do
# Split each element in $a (each line in $1) into an array
	IFS=' ' read -r -a array <<< ${a[$i]}
	### Detection of TE insertion loci (IL)
	# Identification of insertion loci supported by split-mapped reads
	# The output bed file format:
		# Column 1: chromosomal ID of IL
	  # Column 2: insertion coordinate - 5
	  # Column 3: insertion coordinate + 5
	  # Column 4: TE name
	  # Column 5: 0 and 1 indicate the read being mapped to the left and right end of TE, respectively
	  # Column 6: + and - indicate the TE was inserted in the forward and reverse directions in the reference, respectively    
if [[ ${array[2]} = "LEFT" ]]; then
	samtools view -q 1 ${SAMP}.TE_masked.srt.bam ${array[0]}:1-10 | \
awk '{if ($0~/SA:Z:/ && $4<=10) 
	{match($0,/SA:Z:([^,]*),([^,]*),([+-])(,[^,]*,)([^,]*),/,a);
	if (a[4]~/,[0-9]+S[0-9]+M,/ && a[5]>0) 
	    print a[1]"\t"a[2]+$4-1-1"\t"a[2]+$4-1"\t"$3"\t0\t-\t"$1;
	if (a[4]~/,[0-9]+M[0-9]+S,/ && a[5]>0) 
	    {match(a[4],/,([0-9]+)M/,b); print a[1]"\t"a[2]+b[1]-$4-1"\t"a[2]+b[1]-$4"\t"$3"\t0\t+\t"$1};
	if (a[4]~/,[0-9]+S[0-9]+M[0-9]+S,/ && a[5]>0) 
	    {
      match(a[4],/,([0-9]+)S([0-9]+)M([0-9]+)S,/,b); 
			if (b[1]>b[3]) 
			   print a[1]"\t"a[2]+$4-1-1"\t"a[2]+$4-1"\t"$3"\t0\t-\t"$1; 
			else 
			   print a[1]"\t"a[2]+b[2]-$4-1"\t"a[2]+b[2]-$4"\t"$3"\t0\t+\t"$1; 
      }
	}}' | awk -v x=$TE '$1!~x' >> ${SAMP}.SPLIT_SUP.bed
fi
if [[ ${array[2]} = "RIGHT" ]]; then
  RCOOR=${array[1]}
	samtools view -q 1 ${SAMP}.TE_masked.srt.bam ${array[0]}:$((RCOOR-10))-${RCOOR} | \
awk -v x=${RCOOR} -v y=$READ_LEN '{if ($0~/SA:Z:/ && $4>x-y+15) 
	{ match($0,/SA:Z:([^,]*),([^,]*),([+-])(,[^,]*,)([^,]*),/,a); {match($6,/([0-9]+)M/,c);
	if (a[4]~/,[0-9]+S[0-9]+M,/ && a[5]>0 && $6!~/M.*M/ && x-$4-c[1]<10) 
	    print a[1]"\t"a[2]+(x-$4)-c[1]-1"\t"a[2]+(x-$4)-c[1]"\t"$3"\t1\t+\t"$1};
	if (a[4]~/,[0-9]+M[0-9]+S,/ && a[5]>0 && $6!~/M.*M/ && x-$4-c[1]<10) 
	    {match(a[4],/,([0-9]+)M([0-9]+)S,/,b); print a[1]"\t"a[2]+b[1]-x+$4+c[1]-1"\t"a[2]+b[1]-x+$4+c[1]"\t"$3"\t1\t-\t"$1};
	if (a[4]~/,[0-9]+S[0-9]+M[0-9]+S,/ && a[5]>0 && $6!~/M.*M/ && x-$4-c[1]<10) 
	    {match(a[4],/,([0-9]+)S([0-9]+)M([0-9]+)S,/,b); 
			if (b[1]>b[3]) 
			   print a[1]"\t"a[2]+(x-$4)-c[1]-1"\t"a[2]+(x-$4)-c[1]"\t"$3"\t1\t+\t"$1; 
			else 
			   print a[1]"\t"a[2]+b[2]-x+$4+c[1]-1"\t"a[2]+b[2]-x+$4+c[1]"\t"$3"\t1\t-\t"$1;
      }
	}}' | awk -v x=$TE '$1!~x' >> ${SAMP}.SPLIT_SUP.bed
fi
i=$((i+1)) 
done

# Write the informative split-mapped reads supporting the ILs
cut -f 7 ${SAMP}.SPLIT_SUP.bed | sort -u > ${SAMP}_INFOSPLIT_READS.txt

# MQ tag (MAPQ of mate read) is added by samblaster
# Disconcordant read pairs with pair-end reads mapped to reverse (R) strand of chromosomes and TE LEFT end, '-' insertion type
samtools view -h ${SAMP}.disc.readname_srt.bam | samblaster --addMateTags | samtools view -f 48 | grep -v "#" | \
awk -v x=$TE 'match($0, /MQ:i:([0-9]+)/,a) && a[1]>0 && $3~x && $7!~x && $7!~/=/ && $5>0' > ${SAMP}.disc.LEFT_R.IL_R.sam
sed -i "s/\t/|/;s/\t/ /g;s/|/\t/" ${SAMP}.disc.LEFT_R.IL_R.sam # Format multiple column SAM into two columns
Filter_Split_Read.py ${SAMP}_INFOSPLIT_READS.txt ${SAMP}.disc.LEFT_R.IL_R.sam # Filter informative split reads
sed -i "s/ /\t/g" SPLIT_FILT.${SAMP}.disc.LEFT_R.IL_R.sam # Turn two-column SAM into ordinary SAM file

# Disconcordant read pairs with pair-end reads mapped to R strand of TE LEFT end and forward (F) strand of chromosomes, '+' insertion type
samtools view -h ${SAMP}.disc.readname_srt.bam | samblaster --addMateTags | samtools view -f 16 -F 32 | grep -v "#" | \
awk -v x=$TE 'match($0, /MQ:i:([0-9]+)/,a) && a[1]>0 && $3~x && $7!~x && $7!~/=/ && $5>0' > ${SAMP}.disc.LEFT_R.IL_F.sam
sed -i "s/\t/|/;s/\t/ /g;s/|/\t/" ${SAMP}.disc.LEFT_R.IL_F.sam # Format multiple column SAM into two columns
Filter_Split_Read.py ${SAMP}_INFOSPLIT_READS.txt ${SAMP}.disc.LEFT_R.IL_F.sam # Filter informative split reads
sed -i "s/ /\t/g" SPLIT_FILT.${SAMP}.disc.LEFT_R.IL_F.sam # Turn two-column SAM into ordinary SAM file

# Disconcordant read pairs with pair-end reads mapped to F strand of chromosomes and TE Right end, '-' insertion type
samtools view -h ${SAMP}.disc.readname_srt.bam | samblaster --addMateTags | samtools view -F 48 | grep -v "#" | \
awk -v x=$TE 'match($0, /MQ:i:([0-9]+)/,a) && a[1]>0 && $3~x && $7!~x && $7!~/=/ && $5>0' > ${SAMP}.disc.RIGHT_F.IL_F.sam
sed -i "s/\t/|/;s/\t/ /g;s/|/\t/" ${SAMP}.disc.RIGHT_F.IL_F.sam # Format multiple column SAM into two columns
Filter_Split_Read.py ${SAMP}_INFOSPLIT_READS.txt ${SAMP}.disc.RIGHT_F.IL_F.sam # Filter informative split reads
sed -i "s/ /\t/g" SPLIT_FILT.${SAMP}.disc.RIGHT_F.IL_F.sam # Turn two-column SAM into ordinary SAM file

# Disconcordant read pairs with pair-end reads mapped to F strand of TE right end and R strand of chromosomes, '+' insertion type
samtools view -h ${SAMP}.disc.readname_srt.bam | samblaster --addMateTags | samtools view -f 32 -F 16 | grep -v "#" | \
	awk -v x=$TE 'match($0, /MQ:i:([0-9]+)/,a) && a[1]>0 && $3~x && $7!~x && $7!~/=/ && $5>0' > ${SAMP}.disc.RIGHT_F.IL_R.sam
sed -i "s/\t/|/;s/\t/ /g;s/|/\t/" ${SAMP}.disc.RIGHT_F.IL_R.sam # Transfer multiple column SAM into two columns
Filter_Split_Read.py ${SAMP}_INFOSPLIT_READS.txt ${SAMP}.disc.RIGHT_F.IL_R.sam # Filter informative split reads
sed -i "s/ /\t/g" SPLIT_FILT.${SAMP}.disc.RIGHT_F.IL_R.sam # Turn two-column SAM into ordinary SAM file

# Identification of insertion loci supported by disconcordantly mapped reads
  # The output bed file format:
	  # Column 1: chromosomal ID of IL
	  # Column 2: minimum possible coordinate of IL
	  # Column 3: maximum possible coordinate of IL 
	  # Column 4: TE name
	  # Column 5: 0 and 1 indicate the read being mapped to the left and right end of TE, respectively
	  # Column 6: + and - indicate the TE was inserted in the forward and reverse directions in the reference, respectively 
# Disconcordant reads surround Left ends of TEs
awk -v y=${MAX_INS} -v z=${READ_LEN} '{if ($4<y-z) print $7"\t"$8-y+z+int(z/2)"\t"$8+5"\t"$3"\t0\t-\t"$1}' \
	SPLIT_FILT.${SAMP}.disc.LEFT_R.IL_R.sam >> ${SAMP}.DISC_SUP.bed
awk -v y=${MAX_INS} -v z=${READ_LEN} '{match($14, /([0-9]+)M[0-9]*[DSHI]*([0-9]*)M*/, a); if ($4<y-z) print $7"\t"$8+a[1]+a[2]-6"\t"$8+y-z-int(z/2)"\t"$3"\t0\t+\t"$1}' \
	SPLIT_FILT.${SAMP}.disc.LEFT_R.IL_F.sam >> ${SAMP}.DISC_SUP.bed
# Disconcordant reads surround Right ends of TEs 
awk -v y=${MAX_INS} -v z=${READ_LEN} '{match($14, /([0-9]+)M[0-9]*[DSHI]*([0-9]*)M*/, a); print $7"\t"$8+a[1]+a[2]-6"\t"$8+y-z-int(z/2)"\t"$3"\t1\t-\t"$1}' \
	SPLIT_FILT.${SAMP}.disc.RIGHT_F.IL_F.sam >> ${SAMP}.DISC_SUP.bed
awk -v y=${MAX_INS} -v z=${READ_LEN} '{print $7"\t"$8-y+z+int(z/2)"\t"$8+5"\t"$3"\t1\t+\t"$1}' \
	SPLIT_FILT.${SAMP}.disc.RIGHT_F.IL_R.sam >> ${SAMP}.DISC_SUP.bed
```

# # IL results statistics

```bash

# Separate MU1 and MU2 ILs
for SAMP in $( cat ALL1014.SAMP ); do 
awk '$4~/Mu1/' ${SAMP}.SPLIT_SUP.bed > MU1/${SAMP}.MU1.SPLIT_SUP.bed
awk '$4~/Mu2/' ${SAMP}.SPLIT_SUP.bed > MU2/${SAMP}.MU2.SPLIT_SUP.bed 
done

# Put all ILs detected by SPLIT reads in one bed file
cat MU1/*.MU1.SPLIT_SUP.bed > ALL_MU1.SPLIT.bed
cat MU2/*.MU2.SPLIT_SUP.bed > ALL_MU2.SPLIT.bed
# Put all ILs detected by DISC reads in one bed file
cat MU1/*.MU1.DISC_SUP.bed > ALL_MU1.DISC.bed
cat MU2/*.MU2.DISC_SUP.bed > ALL_MU2.DISC.bed
# Put all ILs detected by SPLIT and DISC reads in one bed file
cat ALL_MU1.SPLIT.bed ALL_MU1.DISC.bed > ALL_MU1.SPLIT_DISC.bed
cat ALL_MU2.SPLIT.bed ALL_MU2.DISC.bed > ALL_MU2.SPLIT_DISC.bed

bedtools sort -i ALL_MU1.SPLIT.bed | bedtools merge -d 20 -i - > ALL_MU1.SPLIT_IL.merged.bed
bedtools sort -i ALL_MU2.SPLIT.bed | bedtools merge -d 20 -i - > ALL_MU2.SPLIT_IL.merged.bed

bedtools intersect -wo -a ALL_MU1.SPLIT_IL.merged.bed  -b ALL_MU1.SPLIT.bed > ALL_MU1.SPLIT_IL.SPLIT.intersect
bedtools intersect -wo -a ALL_MU1.SPLIT_IL.merged.bed  -b ALL_MU1.DISC.bed > ALL_MU1.SPLIT_IL.DISC.intersect

bedtools intersect -wo -a ALL_MU2.SPLIT_IL.merged.bed  -b ALL_MU2.SPLIT.bed > ALL_MU2.SPLIT_IL.SPLIT.intersect
bedtools intersect -wo -a ALL_MU2.SPLIT_IL.merged.bed  -b ALL_MU2.DISC.bed > ALL_MU2.SPLIT_IL.DISC.intersect

bedtools intersect -c -a ALL_MU1.SPLIT_IL.merged.bed -b ALL_MU1.SPLIT_DISC.bed > ALL_MU1.SPLIT_IL.SPLIT_DISC.count.tsv
awk '$4>1' ALL_MU1.SPLIT_IL.SPLIT_DISC.count.tsv | cut -f 1,2,3 > ALL_MU1.2SUP.bed

bedtools intersect -c -a ALL_MU2.SPLIT_IL.merged.bed -b ALL_MU2.SPLIT_DISC.bed > ALL_MU2.SPLIT_IL.SPLIT_DISC.count.tsv
awk '$4>1' ALL_MU2.SPLIT_IL.SPLIT_DISC.count.tsv | cut -f 1,2,3 > ALL_MU2.2SUP.bed

# Screen ILs supported by DISC read pairs but not SPLIT read pairs 
bedtools sort -i ALL_MU1.DISC.bed | bedtools intersect -v -wo -a - -b ALL_MU1.SPLIT_IL.merged.bed > ALL_MU1.DISC.no_SPLIT.bed
bedtools sort -i ALL_MU2.DISC.bed | bedtools intersect -v -wo -a - -b ALL_MU2.SPLIT_IL.merged.bed > ALL_MU2.DISC.no_SPLIT.bed
# Merge ILs supported by DISC read pairs but not by SPLIT read pairs, and filter the obtained IL based on the count of supporting read pairs.
bedtools merge -i ALL_MU1.DISC.no_SPLIT.bed | bedtools intersect -wo -a - -b ALL_MU1.DISC.no_SPLIT.bed | awk '$4>1' > ALL_MU1.DISC.no_SPLIT.DISC_intersect.bed
bedtools merge -i ALL_MU2.DISC.no_SPLIT.bed | bedtools intersect -wo -a - -b ALL_MU2.DISC.no_SPLIT.bed | awk '$4>1' > ALL_MU2.DISC.no_SPLIT.DISC_intersect.bed
# Count no_SPLIT IL read support count
bedtools merge -i ALL_MU1.DISC.no_SPLIT.bed | bedtools intersect -c -a - -b ALL_MU1.DISC.no_SPLIT.bed > ALL_MU1.DISC_ONLY_IL.DISC_COUNT.TSV
bedtools merge -i ALL_MU2.DISC.no_SPLIT.bed | bedtools intersect -c -a - -b ALL_MU2.DISC.no_SPLIT.bed > ALL_MU2.DISC_ONLY_IL.DISC_COUNT.TSV
# DISC support only ILs
cat ALL_MU1.DISC_ONLY_IL.DISC_COUNT.TSV | awk '$4>1 {print $1"\t"$2"\t"$3}' > ALL_MU1.DISC_ONLY_IL.bed
cat ALL_MU2.DISC_ONLY_IL.DISC_COUNT.TSV | awk '$4>1 {print $1"\t"$2"\t"$3}' > ALL_MU2.DISC_ONLY_IL.bed
# Count reads in each sample and put all the data into one table file
cp ALL_MU1.DISC_ONLY_IL.bed MU1.DISC_ONLY.244SAMP.READ_COUNT.tsv
for SAMP in $( cat 1124_ALL_SAMP.list ); do
bedtools intersect -c -a ALL_MU1.DISC_ONLY_IL.bed -b MU1/${SAMP}.MU1.DISC_SUP.bed | cut -f 4 > TEMP_COLUMN.txt
paste MU1.DISC_ONLY.244SAMP.READ_COUNT.tsv TEMP_COLUMN.txt > MU1.DISC_ONLY.244SAMP.READ_COUNT.tsv.tmp
mv MU1.DISC_ONLY.244SAMP.READ_COUNT.tsv.tmp MU1.DISC_ONLY.244SAMP.READ_COUNT.tsv
done
cp ALL_MU2.DISC_ONLY_IL.bed MU2.DISC_ONLY.244SAMP.READ_COUNT.tsv
for SAMP in $( cat 1124_ALL_SAMP.list ); do
bedtools intersect -c -a ALL_MU2.DISC_ONLY_IL.bed -b MU2/${SAMP}.MU2.DISC_SUP.bed | cut -f 4 > TEMP_COLUMN.txt
paste MU2.DISC_ONLY.244SAMP.READ_COUNT.tsv TEMP_COLUMN.txt > MU2.DISC_ONLY.244SAMP.READ_COUNT.tsv.tmp
mv MU2.DISC_ONLY.244SAMP.READ_COUNT.tsv.tmp MU2.DISC_ONLY.244SAMP.READ_COUNT.tsv
done

tr '\n' '\t' < 1124_ALL_SAMP.list > header.txt
echo -e "chr\tLeft\tRight\t$( cat header.txt )" > header.txt
cat header.txt MU1.DISC_ONLY.244SAMP.READ_COUNT.tsv > MU1.DISC_ONLY.244SAMP.READ_COUNT.tsv.tmp
mv MU1.DISC_ONLY.244SAMP.READ_COUNT.tsv.tmp MU1.DISC_ONLY.244SAMP.READ_COUNT.tsv
cat header.txt MU2.DISC_ONLY.244SAMP.READ_COUNT.tsv > MU2.DISC_ONLY.244SAMP.READ_COUNT.tsv.tmp
mv MU2.DISC_ONLY.244SAMP.READ_COUNT.tsv.tmp MU2.DISC_ONLY.244SAMP.READ_COUNT.tsv
```

## Identify the MU1 and MU2 members in SWO assemblies

```bash
for DB in GCA_019144185.1 GCA_019144245.1 GCA_019144225.1 GCA_019144155.1 GCA_019144195.1 GCA_019143665.1 GCA_018104345.1 CSIV4 ; do
makeblastdb -in ${DB}.fasta -dbtype nucl -parse_seqids -out $DB
THREADS=8
blastn -query DVS_Mu.fasta -db ${DB} -outfmt 6 -num_threads $THREADS | \
	awk '$3>95 && $4>=500' | cut -f 2,9,10 | \
	awk '$2<$3 {print $1"\t"$2"\t"$3} $2>$3 {print $1"\t"$3"\t"$2}' | \
	bedtools sort -i - | bedtools merge -i - > ${DB}.MU.TEMP.bed

blastn -query DVS_Mu.fasta -db ${DB} -outfmt 6 -num_threads $THREADS | \
	awk '$3>95 && $4>=500 && $9<$10 {print $2"\t"$9"\t"$10"\t"$1"\t"$12"\t+"} $3>95 && $4>=500 && $9>$10 {print $2"\t"$10"\t"$9"\t"$1"\t"$12"\t-"}' > DVS_MU.${DB}.score.bed

bedtools intersect -wo -a ${DB}.MU.TEMP.bed -b DVS_MU.${DB}.score.bed | awk '{print $1":"$2"-"$3"\t"$0}' | sort -rgk 9,9 | awk '!a[$1]++' > ${DB}.MU.BESTHIT.bed

awk '$11>5000 {print $2"\t"$3-500"\t"$3"\n"$2"\t"$4"\t"$4+500}' ${DB}.MU.BESTHIT.bed > ${DB}.MU.bed

bedtools getfasta -fi ${DB}.fasta -fo ${DB}.MU.UP_DOWN_500.fasta -bed ${DB}.MU.bed
blastn -num_threads ${THREADS} -query ${DB}.MU.UP_DOWN_500.fasta -db /zfs/socbd/bwu4/YU/annotation/CK2021 -outfmt 6 | awk '$3>95 && $4>400' > ${DB}.MU.UP_DOWN_500.CK2021.results 
done
```

```bash
### Detect SVs in SWO assemblies
cd /zfs/socbd/bwu4/YU/annotation/FOUR_TE/TE_IN_CITRUS_GENOMES/SWO
ln /zfs/socbd/bwu4/YU/yu_pacbio/SV/MUMMER/T19.asem.fasta
ln /zfs/socbd/bwu4/YU/yu_pacbio/SV/MUMMER/SF.asem.fasta
ln /zfs/socbd/bwu4/YU/yu_pacbio/SV/MUMMER/T78.asem.fasta
for FASTA in Csiv4.chromosome.fa GCA_018104345.1.fasta GCA_019143665.1.fasta GCA_019144155.1.fasta GCA_019144185.1.fasta GCA_019144195.1.fasta GCA_019144225.1.fasta GCA_019144245.1.fasta SF.asem.fasta T78.asem.fasta ; do
minimap2 -cx asm5 -t20 --cs CK2021.60.corrected.fasta $FASTA  > asm.paf  # keeping this file is recommended; --cs required!
sort -k6,6 -k8,8n asm.paf > asm.srt.paf             # sort by reference start coordinate
k8 paftools.js call asm.srt.paf > ${FASTA}.var.txt
done
for FASTA in Csiv4.chromosome.fa GCA_018104345.1.fasta GCA_019143665.1.fasta GCA_019144155.1.fasta GCA_019144185.1.fasta GCA_019144195.1.fasta GCA_019144225.1.fasta GCA_019144245.1.fasta SF.asem.fasta T78.asem.fasta ; do
cat ${FASTA}.var.txt | awk '$5==1 && ($4-$3>=50 || $11-$10>=50) && $4-$3<20000 && $11-$10<20000 { if ($4-$3>=50) {print ">"$2"_"$3"_"$4"\t"$7} else {print ">"$2"_"$3"_"$4"\t"$8}}' >> LARGE_INDEL.tsv
done
awk '!a[$1]++' LARGE_INDEL.tsv | sed "s/\t/\n/g" > LARGE_INDEL.fasta
cd-hit-est -r 1 -g 1 -c 0.95 -i LARGE_INDEL.fasta -o LARGE_INDEL.95_90.clusters.fasta -T 0 -aL 0.90 -M 120000 -d 50 -n 10
cat LARGE_INDEL.95_90.clusters.fasta.clstr | awk '{if ($0~/>Cluster/ && c>3) {print a"\t"b"\t"c;a=$0;b=0;c=0} else if ($0~/*/) {b=$0;c+=1} else {c+=1}}' | sort -n -k7
```

# Check IL in Pacbio CLR reads

```bash
# Generate a new masked reference with separate intact DVS MU clstr sequences
# Put all MU member coordinates in the file DVS_MU.bed
cd /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW
bedtools maskfasta -fi /zfs/socbd/bwu4/YU/annotation/CK2021.60.corrected.fasta -bed DVS_MU.bed -fo CK2021.MU_CLSTR_SEPARATE.fasta
cat /zfs/socbd/bwu4/YU/yu_pacbio/Citrus_sinensis_CK/CORRECTED_DVS_Mu_CLSTR.fasta >> CK2021.MU_CLSTR_SEPARATE.fasta

# To accelerate separating the MU overlapping reads, mapping all raw CLR reads to MU CLSTR sequences: CORRECTED_DVS_Mu_CLSTR.fasta
minimap2 -t 36 -ax map-pb CORRECTED_DVS_Mu_CLSTR.fasta CK.raw.fasta | samtools view -h -F 4 > CK_CLR.MU_CLSTR.sam
samtools fasta CK_CLR.MU_CLSTR.sam > CK_CLR.MU_CLSTR.fasta

# Mapping MU overlapping reads to the masked reference with separate intact DVS MU clstr sequences
minimap2 -t 36 -ax map-pb /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/CK2021.MU_CLSTR_SEPARATE.fasta CK_CLR.MU_CLSTR.fasta | samtools sort -o CK.CK2021_MU_CLSTR.srt.bam

for SAMP in T19 SF T78 ; do
cd /zfs/socbd/bwu4/YU/yu_pacbio/Citrus_sinensis_${SAMP}
minimap2 -t 36 -ax map-pb /zfs/socbd/bwu4/YU/yu_pacbio/Citrus_sinensis_CK/CORRECTED_DVS_Mu_CLSTR.fasta /zfs/socbd/bwu4/YU/yu_pacbio/PACBIO_FASTQ/${SAMP}.fastq.gz | samtools view -h -F 4 > ${SAMP}_CLR.MU_CLSTR.sam
samtools fasta ${SAMP}_CLR.MU_CLSTR.sam > ${SAMP}_CLR.MU_CLSTR.fasta
minimap2 -t 36 -ax map-pb /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/CK2021.MU_CLSTR_SEPARATE.fasta ${SAMP}_CLR.MU_CLSTR.fasta | samtools sort -o ${SAMP}.CK2021_MU_CLSTR.srt.bam
done

samtools view ${SAMP}.CK2021_MU_CLSTR.srt.bam | \
awk '$5>0 && $3~/Mu/ && $0~/SA:Z/ {match($0,/(SA:Z:[^;]+);/,a); if (a[1]!~/Mu/) print $1"\t"$3"\t"$4"\t"a[1]}' > ${SAMP}.temp.tsv
cat ${SAMP}.temp.tsv | awk '{match($1,/\/([0-9]+)_([0-9]+)/,b); match($4,/SA:Z:([^,]+),([^,]+),([^,]+),([^,]+),([^,]+),/,a) 
if (a[4]~/[0-9]+S[0-9A-Z]+[0-9]S/ && a[5]>0) 
	{
	match(a[4],/([0-9]+)S([0-9]+)M[0-9]+I([0-9]+)S/,c); 
	if (c[1]>c[3]) 
	print a[1]"\t"a[2]-100"\t"a[2]+100"\t"$2
	else
	print a[1]"\t"a[2]+c[2]-100"\t"a[2]+c[2]+100"\t"$2
	}
else if (a[4]~/[0-9A-Z]+[0-9]S/ && a[5]>0) 
	{
	match(a[4],/([0-9]+)M[0-9]+I([0-9]+)S/,c);
	print a[1]"\t"a[2]+c[1]-100"\t"a[2]+c[1]+100"\t"$2
	}
else if (a[4]~/([0-9]+)S[0-9A-Z]+/ && a[5]>0)
	{
  print a[1]"\t"a[2]-100"\t"a[2]+100"\t"$2
  } 
}' | awk '$2>=0 {print $0} $2<0 {print $1"\t0\t"$3"\t"$4}' | bedtools sort -i - > ${SAMP}.MU_SORT.bed
bedtools merge -i ${SAMP}.MU_SORT.bed > ${SAMP}.MU_MERGE.bed
bedtools intersect -c -a ${SAMP}.MU_MERGE.bed -b ${SAMP}.MU_SORT.bed | awk '$4>3' > ${SAMP}.MU_SORT.4ZMW.bed
bedtools intersect -v -a ${SAMP}.MU_SORT.4ZMW.bed -b ${SAMP}.SPLIT_SUP.CLR_NGS.COUNT.tsv > ${SAMP}.temp.txt
bedtools intersect -v -a ${SAMP}.temp.txt -b /zfs/socbd/bwu4/YU/annotation/FOUR_TE/TE_IN_CITRUS_GENOMES/CORRECTED_DVS_4TE.bed > ${SAMP}.CLR_SPECIFIC.txt

cat /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/${SAMP}.SPLIT_SUP.bed /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/${SAMP}.DISC_SUP.bed > /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/${SAMP}.MERGE.bed
bedtools sort -i /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/${SAMP}.SPLIT_SUP.bed | \
bedtools merge -d 15 -i - | bedtools intersect -c -a - -b ${SAMP}.MU_SORT.bed | \
bedtools intersect -c -a - -b /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/${SAMP}.MERGE.bed > ${SAMP}.SPLIT_SUP.CLR_NGS.COUNT.tsv
```

# Transcripts from MU2 in RNAseq data

```bash
makeblastdb -in /zfs/socbd/bwu4/YU/annotation/CHRA/chrA_mikado.tr.fasta -dbtype nucl -parse_seqids -out CHRA_MIKADO_TR
makeblastdb -in /zfs/socbd/bwu4/YU/annotation/CHRB/chrB_mikado.tr.fasta -dbtype nucl -parse_seqids -out CHRB_MIKADO_TR

blastn -query CORRECTED_DVS_Mu_CLSTR.fasta -db CHRA_MIKADO_TR -num_threads 8 -outfmt 6 \
| awk '$3>95 && $4>1000 && $1~/Mu2/' | cut -f 2 | sort -u | \
cdbyank  /zfs/socbd/bwu4/YU/annotation/CHRA/chrA_mikado.tr.fasta.cidx > RNAseq_asembly.MU2.fasta

blastn -query CORRECTED_DVS_Mu_CLSTR.fasta -db CHRB_MIKADO_TR -num_threads 8 -outfmt 6 \
| awk '$3>95 && $8>6100 && $1=="DVS_Mu2_1"' | cut -f 2 | sort -u > chrB.MU2.ID.txt

/zfs/socbd/bwu4/YU/annotation/CHRA/MIKDAO_CHRA/5-mikado-permissive/pick/permissive/

blastn -query CORRECTED_DVS_Mu_CLSTR.fasta -db chrB_PERMISSIVE -num_threads 8 -outfmt 6 \
| awk '$3>95 && $8>6100 && $1=="DVS_Mu2_1"' | cut -f 2 | sort -u > chrB.MU2.ID.txt

gffread /zfs/socbd/bwu4/YU/annotation/CHRA/MIKDAO_CHRA/5-mikado-lenient/pick/permissive/mikado-permissive.loci.gff3 -g /zfs/socbd/bwu4/YU/annotation/CHRA/chrA.fasta -w chrA_mikado.permissive.tr.fasta
gffread /zfs/socbd/bwu4/YU/annotation/CHRA/MIKDAO_CHRA/5-mikado-lenient/pick/lenient/mikado-lenient.loci.gff3 -g /zfs/socbd/bwu4/YU/annotation/CHRA/chrA.fasta -w chrA_mikado.lenient.tr.fasta

makeblastdb -in chrA_mikado.permissive.tr.fasta -dbtype nucl -parse_seqids -out chrA_PERMISSIVE
makeblastdb -in chrA_mikado.lenient.tr.fasta -dbtype nucl -parse_seqids -out chrA_LENIENT
blastn -query CORRECTED_DVS_Mu_CLSTR.fasta -db chrA_LENIENT -num_threads 8 -outfmt 6 | awk '$3>95 && $8>6100 && $1=="DVS_Mu2_1" && $9>$10'
blastn -query CORRECTED_DVS_Mu_CLSTR.fasta -db chrA_LENIENT -num_threads 8 -outfmt 6 \
| awk '$3>95 && $8>6100 && $1=="DVS_Mu2_1" && $9>$10' | cut -f 2 | sort -u > chrA.MU2.ID.txt

gffread /zfs/socbd/bwu4/YU/annotation/CHRB/MIKADO/5-mikado-lenient/pick/lenient/mikado-lenient.loci.gff3 -g /zfs/socbd/bwu4/YU/annotation/CHRB/chrB.fasta -w chrB_mikado.lenient.tr.fasta
makeblastdb -in chrB_mikado.lenient.tr.fasta -dbtype nucl -parse_seqids -out chrB_LENIENT
blastn -query CORRECTED_DVS_Mu_CLSTR.fasta -db chrB_LENIENT -num_threads 8 -outfmt 6 \
| awk '$3>95 && $8>6100 && $1=="DVS_Mu2_1" && $9>$10' | cut -f 2 | sort -u > chrB.MU2.ID.txt

cdbfasta chrA_mikado.lenient.tr.fasta
cdbfasta chrB_mikado.lenient.tr.fasta
cat chrA.MU2.ID.txt | cdbyank chrA_mikado.lenient.tr.fasta.cidx - > chrA.MU2.LENIENT.fasta
cat chrB.MU2.ID.txt | cdbyank chrB_mikado.lenient.tr.fasta.cidx - > chrB.MU2.LENIENT.fasta
```

# Comparison between SWO and PTR Mu sequences

```bash
sed -i 's/(+)//g;s/(-)//g' CORRECTED_DVS_Mu_CLSTR.fasta
cdbfasta CORRECTED_DVS_Mu_CLSTR.fasta
for MU in $( grep ">" CORRECTED_DVS_Mu_CLSTR.fasta | sed "s/>//g" ); do cdbyank CORRECTED_DVS_Mu_CLSTR.fasta.cidx -a $MU > ${MU}.fasta ; done
sed -i 's/(-)//g;s/(+)//g' PTR_FOUR_TE.fasta
cdbfasta PTR_FOUR_TE.fasta
for MU in $( grep ">" PTR_FOUR_TE.fasta | sed "s/>//g" ); do cdbyank PTR_FOUR_TE.fasta.cidx -a $MU > ${MU}.fasta ; done
mkdir DVS_PTR
mv *_Mu[12]_*.fasta DVS_PTR
cd DVS_PTR
for DVS in DVS_Mu1*.fasta ; do for PTR in PTR_Mu1*.fasta ; do lalign36 $DVS $PTR 5 > ${DVS}_${PTR}.align.txt ; done ; done

```

## PAIRWISE_ALIGNMENT_DISTANCE.py Pairwise sequence identity calculation

```python
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import sys
import os
aln = AlignIO.read(sys.argv[2],"fasta")
calculator = DistanceCalculator(sys.argv[1])
dm = calculator.get_distance(aln)
print(dm.names[0],dm.names[0],dm[0][1],sep="\t")
```

## Calculate pairwise identity among transposons and genes

```bash
ls -1 PTR_Mu1_?.fasta > PTR_Mu1.list
ls -1 PTR_Mu2_?.fasta > PTR_Mu2.list
ls -1 DVS_Mu2_?.fasta DVS_Mu2_??.fasta > DVS_Mu2.list
ls -1 DVS_Mu1_?.fasta DVS_Mu1_??.fasta > DVS_Mu1.list
for ONE in $( cat DVS_Mu1.list ) ; do
for TWO in $( cat PTR_Mu1.list ) ; do
cat $ONE $TWO > ${ONE%.fasta}_$TWO
muscle -in ${ONE%.fasta}_$TWO -out ${ONE%.fasta}_${TWO%fasta}aligned.fasta 
trimal -nogaps -in ${ONE%.fasta}_${TWO%fasta}aligned.fasta > ${ONE%.fasta}_${TWO%fasta}aligned.nogaps.fasta
./PAIRWISE_ALIGNMENT_DISTANCE.py identity ${ONE%.fasta}_${TWO%fasta}aligned.nogaps.fasta >> PTR_DVS_MU_DIST.txt
done
done
```

## Genetic distance control group calculation

```bash
cd /zfs/socbd/bwu4/YU/annotation/ORTHOFINDER/ORTHO07012021/Results_Jul01/Single_Copy_Orthologue_Sequences
for FA in *.fa ; do grep ">" $FA | sed "s/>//g" | \
awk '{ if ($0~/DSWO[0-9]A/ || $0~/Ptrif/) print $0}' \
> ${FA%fa}DVSA_PTR.txt; done
for FA in *.fa ; do grep ">" $FA | sed "s/>//g" | \
awk '{ if ($0~/DSWO[0-9]B/ || $0~/Ptrif/) print $0}' \
> ${FA%fa}DVSB_PTR.txt; done
cd /zfs/socbd/bwu4/REF/DVS
awk '$3~/transcript/ && match($9,/ID=([^;]*);/,a) {print $1"\t"$4"\t"$5"\t"a[1]"\t.\t"$7 }' \
CK2021.FINAL.0403.gff3 > CK2021.tr_genomic.bed
bedtools getfasta -s -fi CK2021.60.corrected.fasta -fo \
CK2021.tr_genomic.fasta -nameOnly -bed CK2021.tr_genomic.bed
mv CK2021.tr_genomic.fasta /zfs/socbd/bwu4/YU/annotation/FOUR_TE/TE_IN_CITRUS_GENOMES/DVS_PTR/CONTROL
awk '$3~/mRNA/ && match($9,/ID=([^;]*);/,a) {print $1"\t"$4"\t"$5"\t"a[1]"\t.\t"$7 }' \
Ptrifoliata_565_v1.3.1.gene.gff3 > PTR.tr_genomic.bed
bedtools getfasta -s -fi Ptrifoliata_565_v1.3.fa -fo PTR.tr_genomic.fasta \
-nameOnly -bed PTR.tr_genomic.bed
mv Ptrifoliata_565_v1.3.fa /zfs/socbd/bwu4/YU/annotation/FOUR_TE/TE_IN_CITRUS_GENOMES/DVS_PTR/CONTROL
cat CK2021.tr_genomic.fasta Ptrifoliata_565_v1.3.fa > DVS_PTR.REF.fasta
sed -i 's/(+)//g;s/(-)//g' DVS_PTR.REF.fasta
cdbfasta DVS_PTR.REF.fasta
ls -1 /zfs/socbd/bwu4/YU/annotation/ORTHOFINDER/ORTHO07012021/Results_Jul01/Single_Copy_Orthologue_Sequences/*.DVSA_PTR.txt | \
head -n 500 > 500OG.DVSA_PTR.txt
ls -1 /zfs/socbd/bwu4/YU/annotation/ORTHOFINDER/ORTHO07012021/Results_Jul01/Single_Copy_Orthologue_Sequences/*.DVSB_PTR.txt | \
head -n 5000 | tail -n 500 > 500OG.DVSB_PTR.txt

for TXT in $( cat 500OG.DVSA_PTR.txt ); do cat $TXT | \
cdbyank DVS_PTR.REF.fasta.cidx - > ${TXT%txt}fasta ; done

for TXT in $( cat 500OG.DVSB_PTR.txt ); do cat $TXT | \
cdbyank DVS_PTR.REF.fasta.cidx - > ${TXT%txt}fasta ; done
mv /zfs/socbd/bwu4/YU/annotation/ORTHOFINDER/ORTHO07012021/Results_Jul01/Single_Copy_Orthologue_Sequences/OG*.DVS*.fasta ./
```

```bash
for OG in OG*.DVS*.fasta ; do
muscle -in $OG -out ${OG%fasta}aligned.fasta
trimal -nogaps -in ${OG%fasta}aligned.fasta > ${OG%fasta}aligned.nogaps.fasta
../PAIRWISE_ALIGNMENT_DISTANCE.py identity  ${OG%fasta}aligned.nogaps.fasta >> CONTROL_DIST.txt
done

awk '{ if ( $4-$3>=4500 && $12>0 ) print $1"\t"$3"\t"$4"\tDVS"NR"\t.\t+\t"$6"\t"$8"\t"$9"\tPTR"NR"\t.\t"$5}' \
PTR.DVS_INTER_GENE.5000.paf > DVS_PTR.INTER_GENE.OVER4500.tsv
```

```bash
# shuf -n 1000 DVS_PTR.INTER_GENE.OVER4500.tsv > DVS_PTR.INTER_GENE.OVER4500.RAN1000.tsv
cut -f 7,8,9,10,11,12 DVS_PTR.INTER_GENE.OVER4500.tsv > PTR.OVER4500.bed

awk '{ if ( $12>0 && $4-$3>500 ) print $6"\t"$8"\t"$9 }' PTR.DVS_INTER_GENE.5000.paf | \
bedtools sort -i - > PTR_COVERED.bed

cut -f 7,8,9,10,11,12 DVS_PTR.INTER_GENE.OVER4500.tsv > PTR.OVER4500.bed

bedtools intersect -wa -c -a PTR.OVER4500.bed -b PTR_COVERED.bed | \
awk '$7==1' | cut -f 4 > PTR.DVS_UNIQ.ID.txt

#/usr/bin/env python3
import pandas as pd
df = pd.read_table('DVS_PTR.INTER_GENE.OVER4500.tsv', header=None)
df1 = pd.read_table('PTR.DVS_UNIQ.ID.txt', header=None)
df.loc[df[10].isin(df1[0])].to_csv('ORG.DVS_PTR.1VS1.tsv', index=False, header=False)
### python code ends

sed -i "s/\tscaffold/\nscaffold/g" ORG.DVS_PTR.1VS1.tsv

split -l 2 ORG.DVS_PTR.1VS1.tsv DVS_PTR.1VS1.

for FILE in DVS_PTR.1VS1.* ; do 
bedtools getfasta -s -nameOnly -fi CK2021.INTER_GENE.5000.PTR.fasta -fo ${FILE}.fas -bed $FILE
done

cd /zfs/socbd/bwu4/YU/annotation/FOUR_TE/TE_IN_CITRUS_GENOMES/DVS_PTR/CONTROL
mv /zfs/socbd/bwu4/REF/PTR/DVS_PTR.1VS1.*.fas ./

for OG in DVS_PTR.1VS1.*.fas ; do
muscle -in $OG -out ${OG%fas}aligned.fasta
trimal -nogaps -in ${OG%fas}aligned.fasta > ${OG%fas}aligned.nogaps.fasta
../PAIRWISE_ALIGNMENT_DISTANCE.py identity  ${OG%fas}aligned.nogaps.fasta >> CONTROL_DIST.txt
done

arrVar=()
for DVS1 in $( cat PTR_Mu2.list ); do
arrVar=( ${arrVar[@]} $DVS1 )
for DVS2 in $( cat PTR_Mu2.list ); do
if [[ " ${arrVar[*]} " =~ " $DVS2 " ]] ; then
continue
else
cat $DVS1 $DVS2 > ${DVS1%.fasta}_${DVS2}
muscle -in ${DVS1%.fasta}_${DVS2} -out ${DVS1%.fasta}_${DVS2%fasta}aligned.fasta
trimal -nogaps -in ${DVS1%.fasta}_${DVS2%fasta}aligned.fasta > ${DVS1%.fasta}_${DVS2%fasta}aligned.nogaps.fasta
./PAIRWISE_ALIGNMENT_DISTANCE.py identity ${DVS1%.fasta}_${DVS2%fasta}aligned.nogaps.fasta >> PTR_DVS_MU_DIST.txt
fi
done
done

arrVar=()
for TE1 in $( cat DVS_PTR_Mu2.list ); do
arrVar=( ${arrVar[@]} $TE1 )
for TE2 in $( cat DVS_PTR_Mu2.list ); do
if [[ " ${arrVar[*]} " =~ " $TE2 " ]] ; then
continue
else
needle -gapopen 10.0 -gapextend 0.5  -brief -outfile ${TE1%.fasta}_${TE2%.fasta}.needle $TE1 $TE2
a=$( grep "# 1" ${TE1%.fasta}_${TE2%.fasta}.needle | awk -F " " '{print $3}' )
b=$( grep "# 2" ${TE1%.fasta}_${TE2%.fasta}.needle | awk -F " " '{print $3}' )
c=$( grep Identity ${TE1%.fasta}_${TE2%.fasta}.needle | tr -s ' ' | awk -F " " '{match($3,/([0-9]+)\/([0-9]+)/,a); print a[1]}' )
d=$( grep Gaps ${TE1%.fasta}_${TE2%.fasta}.needle | tr -s ' ' | awk -F " " '{match($3,/([0-9]+)\/([0-9]+)/,a); print a[1]}' )
e=$( grep Identity ${TE1%.fasta}_${TE2%.fasta}.needle | tr -s ' ' | awk -F " " '{match($3,/([0-9]+)\/([0-9]+)/,a); print a[2]}' )
echo -e $a"\t"$b"\t"$c"\t"$d"\t"$e >> NEEDLE_MU2_DIST.tsv
fi
done
done

for FAS in DVS_PTR.1VS1.*.fas ; do
split -l 2 -d $FAS ${FAS%.fas}
needle -gapopen 10.0 -gapextend 0.5 -brief -outfile ${FAS%.fas}.needle ${FAS%.fas}00 ${FAS%.fas}01
a=$( grep "# 1" ${FAS%.fas}.needle | awk -F " " '{print $3}' )
b=$( grep "# 2" ${FAS%.fas}.needle | awk -F " " '{print $3}' )
c=$( grep Identity ${FAS%.fas}.needle | tr -s ' ' | awk -F " " '{match($3,/([0-9]+)\/([0-9]+)/,a); print a[1]}' )
d=$( grep Gaps ${FAS%.fas}.needle | tr -s ' ' | awk -F " " '{match($3,/([0-9]+)\/([0-9]+)/,a); print a[1]}' )
e=$( grep Identity ${FAS%.fas}.needle | tr -s ' ' | awk -F " " '{match($3,/([0-9]+)\/([0-9]+)/,a); print a[2]}' )
echo -e $a"\t"$b"\t"$c"\t"$d"\t"$e >> NEEDLE_CONTROL_DIST.tsv
done

TE1="PTR_Mu1_3.fasta"
TE2="PTR_Mu1_5.fasta"
needle -gapopen 10.0 -gapextend 0.5  -brief -outfile ${TE1%.fasta}_${TE2%.fasta}.needle $TE1 $TE2
a=$( grep "# 1" ${TE1%.fasta}_${TE2%.fasta}.needle | awk -F " " '{print $3}' )
b=$( grep "# 2" ${TE1%.fasta}_${TE2%.fasta}.needle | awk -F " " '{print $3}' )
c=$( grep Identity ${TE1%.fasta}_${TE2%.fasta}.needle | tr -s ' ' | awk -F " " '{match($3,/([0-9]+)\/([0-9]+)/,a); print a[1]}' )
d=$( grep Gaps ${TE1%.fasta}_${TE2%.fasta}.needle | tr -s ' ' | awk -F " " '{match($3,/([0-9]+)\/([0-9]+)/,a); print a[1]}' )
e=$( grep Identity ${TE1%.fasta}_${TE2%.fasta}.needle | tr -s ' ' | awk -F " " '{match($3,/([0-9]+)\/([0-9]+)/,a); print a[2]}' )
echo -e $a"\t"$b"\t"$c"\t"$d"\t"$e
```

# Kmer statistics in the reference genomes

```bash
cd /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/KMER_ANALYSIS
for KMER in {47..64}; do 
jellyfish count -C -m $KMER -s 600M -t 24 ../CK2021.09282021.MASKED.fasta 
mv mer_counts.jf mer${KMER}_counts.jf 
done

for KMER in {65..75}; do 
jellyfish count -C -m $KMER -s 600M -t 24 ../CK2021.09282021.MASKED.fasta 
mv mer_counts.jf mer${KMER}_counts.jf 
done

for KMER in {86..91}; do 
jellyfish count -C -m $KMER -s 600M -t 16 ../CK2021.09282021.MASKED.fasta 
mv mer_counts.jf mer${KMER}_counts.jf 
done

KMER=65
jellyfish histo mer${KMER}_counts.jf | sed "s/ /\t/g" > mer${KMER}.count
awk 'BEGIN {a=0} NR==1 {b=$2} {a+=$1*$2} END {print b"\t"a"\t"b/a}' mer${KMER}.count

for KMER in {20..100}; do
if test ! -f mer${KMER}.count; then
jellyfish histo mer${KMER}_counts.jf | sed "s/ /\t/g" > mer${KMER}.count
fi
echo mer${KMER}.count $( awk 'BEGIN {a=0} NR==1 {b=$2} {a+=$1*$2} END {print b"\t"a"\t"b/a}' mer${KMER}.count) >> UNIQ_KMER_STA.txt
done
```

# SWO assembly IL SUL statistics

### IL_SUL_ANALYSIS.sh

```bash
#!/bin/bash

for LEN in {20..100}; do
b=$( sed -n ${1}p ${2} | awk -v x=${LEN} '$3-$2>50 {print $1"\t"$2-x-1"\t"$2-1} $3-$2<50 {print $1"\t"$3-x"\t"$3}' | \
bedtools getfasta -fi /zfs/socbd/bwu4/YU/annotation/CK2021.60.corrected.fasta -bed - | tail -n 1)
a=$( jellyfish query mer${LEN}_counts.jf $b | cut -d " " -f 2 )
echo $a
if test $a -eq 1 ; then
sed -n ${1}p ${2} | awk -v x=${LEN} '{print $1":"$2"-"$3"\t"x"\tLEFT"}' >> ${3}
break
fi
if test $a -eq 0 ; then
sed -n ${1}p ${2} | awk -v x=${LEN} '{print $1":"$2"-"$3"\t"x"_0\tLEFT"}' >> ${3}
break
fi
if test $LEN -eq 100 ; then
sed -n ${1}p ${2} | awk -v x=${LEN} '{print $1":"$2"-"$3"\t"x"_2\tLEFT"}' >> ${3}
fi
done
for LEN in {20..100}; do
b=$( sed -n ${1}p ${2} | awk -v x=${LEN} '$3-$2>50 {print $1"\t"$3"\t"$3+x} $3-$2<50 {print $1"\t"$2-1"\t"$2+x-1}' | \
bedtools getfasta -fi /zfs/socbd/bwu4/YU/annotation/CK2021.60.corrected.fasta -bed - | tail -n 1 )
a=$( jellyfish query mer${LEN}_counts.jf $b | cut -d " " -f 2 )
echo $a
if test $a -eq 1 ; then
sed -n ${1}p ${2} | awk -v x=${LEN} '{print $1":"$2"-"$3"\t"x"\tRIGHT"}' >> ${3}
break
fi
if test $a -eq 0 ; then
sed -n ${1}p ${2} | awk -v x=${LEN} '{print $1":"$2"-"$3"\t"x"_0\tRIGHT"}' >> ${3}
break
fi
if test $LEN -eq 100 ; then
sed -n ${1}p ${2} | awk -v x=${LEN} '{print $1":"$2"-"$3"\t"x"_2\tRIGHT"}' >> ${3}
fi
done
```

```bash
for N in {1..126}; do
./IL_SUL_ANALYSIS.sh $N SWO_ASSEM_IL.txt IL_SUL.tsv
done

echo -e "chr5B\t48342312\t48342412" |  bedtools getfasta -fi /zfs/socbd/bwu4/YU/annotation/CK2021.60.corrected.fasta -bed -
```

### MU END SUL statistics

```bash
for i in {1..22}; do
for LEN in {20..100}; do
b=$( sed -n ${i}p MU_END.txt | awk -v x=${LEN} '$2==1 {print $1"\t0\t"x} $2>1 {print $1"\t"$2-x"\t"$2}' | \
bedtools getfasta -fi ../CK2021.09282021.MASKED.fasta -bed - | tail -n 1)
a=$( jellyfish query mer${LEN}_counts.jf $b | cut -d " " -f 2 )
echo $a
if test $a -eq 1 ; then
sed -n ${i}p MU_END.txt | awk -v x=${LEN} '{print $1"\t"$2"\t"x}' >> MU_END_SUL.tsv
break
fi
if test $a -eq 0 ; then
sed -n ${i}p MU_END.txt | awk '{print $1"\t"$2"\t0"}' >> MU_END_SUL.tsv
fi
if test $LEN -eq 100 ; then
sed -n ${i}p MU_END.txt | awk '{print $1"\t"$2"\t101"}' >> MU_END_SUL.tsv
fi
done
done
```

# Fix the overlapping between multiple SPLIT IL with a single DISC IL

```bash
cd /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW
cat *.DISC_SUP.bed | grep Mu2 > ALL.MU2_DISC.bed
cat *.DISC_SUP.bed | grep Mu1 > ALL.MU1_DISC.bed

bedtools intersect -wo -a ALL.MU2_DISC.UNIQ.bed -b ALL_MU2.2SUP.bed | cut -f 1,2,3 \
| sort | uniq -d | bedtools sort -i - | bedtools merge -i - > MU2.DISC.2SPLIT_IL.bed
bedtools intersect -wo -a ALL.MU1_DISC.UNIQ.bed -b ALL_MU1.2SUP.bed | cut -f 1,2,3 \
| sort | uniq -d | bedtools sort -i - | bedtools merge -i - > MU1.DISC.2SPLIT_IL.bed

bedtools intersect -wo -a MU1.DISC.2SPLIT_IL.bed -b ALL_MU1.2SUP.bed > MU1.DISC.2SPLIT.IL.intersect
bedtools intersect -wo -a MU2.DISC.2SPLIT_IL.bed -b ALL_MU2.2SUP.bed > MU2.DISC.2SPLIT.IL.intersect

cut -f 4,5,6 MU1.DISC.2SPLIT.IL.intersect > MU1.DISC.2SPLIT.IL.bed
cut -f 4,5,6 MU2.DISC.2SPLIT.IL.intersect > MU2.DISC.2SPLIT.IL.bed
```

# Mutant-type cell ratio analysis

```bash
# Add masked intact and non-intact MULE regions in DVS_MU.bed
# Sort the DVS_MU.bed file
bedtools sort -i DVS_MU.bed > DVS_MU.including_nonintact.sorted.bed

# Extend 20 bp at the two terminals of DVS MULEs
awk '{print $1"\t"$2-20"\t"$3+20}' DVS_MU.including_nonintact.sorted.bed > DVS_MU.including_nonintact.sorted.EXT20.bed

# Remove DVS MU ILs from all the MU2 and MU1 ILs 
bedtools intersect -v -wa -a ALL_MU2.2SUP.CORRECTED.bed -b DVS_MU.including_nonintact.sorted.EXT20.bed > ALL_MU2.2SUP.nonDVS.bed
bedtools intersect -v -wa -a ALL_MU1.2SUP.CORRECTED.bed -b DVS_MU.including_nonintact.sorted.EXT20.bed > ALL_MU1.2SUP.nonDVS.bed

# Output surrounding 1000 bp sequences of filtered MULE ILs and separate those on DVS_A and DVS_B
awk '$1~/A/ {print $1"\t"$2-500"\t"$3+500}' ALL_MU1.2SUP.nonDVS.bed > ALL_MU1.2SUP.nonDVS.chrA.SURR500.bed

awk '$1~/A/ {print $1"\t"$2-500"\t"$3+500}' ALL_MU2.2SUP.nonDVS.bed > ALL_MU2.2SUP.nonDVS.chrA.SURR500.bed

# Map the surrounding sequences of DVS_A ILs to unmasked DVS_B to check if their allelic regions are 
# masked in the modified DVS reference genome
bedtools getfasta -fi /zfs/socbd/bwu4/YU/annotation/CK2021.60.corrected.fasta \
-fo ALL_MU1.2SUP.nonDVS.chrA.SURR500.fasta -bed ALL_MU1.2SUP.nonDVS.chrA.SURR500.bed
bedtools getfasta -fi /zfs/socbd/bwu4/YU/annotation/CK2021.60.corrected.fasta \
-fo ALL_MU2.2SUP.nonDVS.chrA.SURR500.fasta -bed ALL_MU2.2SUP.nonDVS.chrA.SURR500.bed
minimap2 -x asm20 -t 8 CK2021.60.corrected.B.fasta ALL_MU1.2SUP.nonDVS.chrA.SURR500.fasta > MU1_chrAIL_SURR500.paf
minimap2 -x asm20 -t 8 CK2021.60.corrected.B.fasta ALL_MU2.2SUP.nonDVS.chrA.SURR500.fasta > MU2_chrAIL_SURR500.paf
cat MU1_chrAIL_SURR500.paf | awk '$12>0 && $4-$3 > 200 {print $6"\t"$8"\t"$9"\t"$1}' > MU1_chrAIL_SURR500_DVSB.bed
cat MU2_chrAIL_SURR500.paf | awk '$12>0 && $4-$3 > 200 {print $6"\t"$8"\t"$9"\t"$1}' > MU2_chrAIL_SURR500_DVSB.bed
bedtools intersect -wo -a MU1_chrAIL_SURR500_DVSB.bed -b CK2021.B.TO_MASK.bed > MU1_chrAIL_SURR500_DVSB_MASKED.intersect.txt
bedtools intersect -wo -a MU2_chrAIL_SURR500_DVSB.bed -b CK2021.B.TO_MASK.bed > MU2_chrAIL_SURR500_DVSB_MASKED.intersect.txt

samtools view -q 1 -f 66 -@ 8 BAM/SRR10150561.TE_masked.srt.bam chr9A:12773268-12773278 | grep -v Mu | \
awk '( $4<=12773268-50 && $8+150>12773278+50) || ($8<=12773268-50 && $4+150>=12773278+50)'

# Prepare reference including DVSA and the modified Mu sequences
cdbfasta CK2021.09282021.MASKED.fasta
grep ">" CK2021.09282021.MASKED.fasta | grep -v B | sed "s/>//g" > DVSA_MU.12022021.gi
cat DVSA_MU.12022021.gi | cdbyank CK2021.09282021.MASKED.fasta.cidx - > DVSA_MU.12022021.fasta
bwa index DVSA_MU.12022021.fasta

for MU in MU1 MU2 ; do
for SAMP in $( cat 1124_ALL_SAMP.list ); do
./FETCH_SAMP_ILS.py ${MU} $SAMP
cat ${SAMP}_${MU}_MASKED.bed | awk '$2>=1000 {print $1"\t"$2-1000"\t"$3+1000} $2<1000 {print $1"\t0\t"$3+1000}' | bedtools sort -i - > ${SAMP}_${MU}_MASK.region
samtools view -h -@ 8 -L ${SAMP}_${MU}_MASK.region BAM/CK.TE_masked.srt.bam |  samtools fastq -1 ${SAMP}_${MU}.TMP_MASK_1.fq -2 ${SAMP}_${MU}.TMP_MASK_2.fq
grep "@" ${SAMP}_${MU}.TMP_MASK_1.fq | sort > ${SAMP}_${MU}.TMP_MASK_1.id
grep "@" ${SAMP}_${MU}.TMP_MASK_2.fq | sort > ${SAMP}_${MU}.TMP_MASK_2.id
comm -1 -2 ${SAMP}_${MU}.TMP_MASK_1.id ${SAMP}_${MU}.TMP_MASK_2.id > ${SAMP}_${MU}.TMP_MASK_12.id
sed -i 's/@//g' ${SAMP}_${MU}.TMP_MASK_12.id
seqtk subseq ${SAMP}_${MU}.TMP_MASK_1.fq ${SAMP}_${MU}.TMP_MASK_12.id | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${MU}.PAIRED.TMP_MASK_1.fq
seqtk subseq ${SAMP}_${MU}.TMP_MASK_2.fq ${SAMP}_${MU}.TMP_MASK_12.id | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${MU}.PAIRED.TMP_MASK_2.fq
bwa mem -t 7 DVSA_MU.12022021.fasta ${SAMP}_${MU}.PAIRED.TMP_MASK_1.fq ${SAMP}_${MU}.PAIRED.TMP_MASK_2.fq | samtools sort -o ${SAMP}_${MU}.PAIRED.TMP_MASK.bam
samtools index ${SAMP}_${MU}.PAIRED.TMP_MASK.bam
rm -f ${SAMP}_${MU}.MASK_W_count.tsv
for IL in $( cat ${SAMP}_${MU}_IL.txt ); do
if [[ ${IL} == *"A"* ]]; then
W=$( samtools view -q 1 -f 66 ${SAMP}_${MU}.PAIRED.TMP_MASK.bam $IL | grep -v Mu | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} ( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${MU}.MASK_W_count.tsv
else
echo -e ${IL}"\t0" >> ${SAMP}_${MU}.MASK_W_count.tsv
fi
done
touch ${SAMP}_${MU}.MASK_W_count.tsv
done
done

# Put all corrected ILs and the genotyeps of the 244 samples in MU1_GENOTYPE.tsv and MU2_GENOTYPE.tsv
for MU in MU1 MU2 ; do
for SAMP in $( cat 1124_ALL_SAMP.list ); do
rm -f ${SAMP}_${MU}.W_count.tsv
for IL in $( cat ${SAMP}_${MU}_IL.txt ); do
ILEXT=$( echo $IL | awk '{match($0,/(chr[1-9AB]+):([0-9]+)-([0-9]+)$/,a); print a[1]":"a[2]-1001"-"a[3]+1000}' )
W=$( samtools view -q 1 -f 66 BAM/${SAMP}.TE_masked.srt.bam $ILEXT | grep -v Mu | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} 
( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${MU}.W_count.tsv
done
touch ${SAMP}_${MU}.W_count.tsv
done
done
```

 # Content of FETCH_SAMP_ILS.py

```python
#!/usr/bin/env python3

import pandas as pd
import sys

MU = sys.argv[1] # Mu name
SAMP = sys.argv[2] # Sample name in the genotype table
GTDF = pd.read_table(f'{MU}_GENOTYPE.tsv', header=0)
GTDF.loc[GTDF[SAMP]==1,'name'].to_csv(f'{SAMP}_{MU}_IL.txt',header=False,index=False)
ALL_BED = pd.read_table(f'{MU}_chrAIL_SURR500_DVSB.bed',names=['chr','Left','Right','surr'])
ID_CONVER = pd.read_table(f'{MU}_IL_ILSURROUNDING.tsv',names=['name','surr'])
ALL_BED = ALL_BED.merge(ID_CONVER,on='surr')
ALL_BED.loc[ALL_BED['name'].isin(GTDF.loc[GTDF[SAMP]==1,'name']),['chr','Left','Right']].to_csv(f'{SAMP}_{MU}_MASKED.bed',sep="\t",header=False,index=False)
ALL_BED.loc[ALL_BED['name'].isin(GTDF.loc[GTDF[SAMP]==1,'name']),['surr']].to_csv(f'{SAMP}_{MU}_SURR.bed',header=False,index=False)

```

 # Merge all ${SAMP}_${MU}.W_count.tsv and ${SAMP}_${MU}.MASK_W_count.tsv

```python
#!/usr/bin/env python3

import pandas as pd

MU1_W_COUNT = pd.read_table('MU1_GENOTYPE.tsv', header=0)
MU1_MASKW_COUNT = pd.read_table('MU1_GENOTYPE.tsv', header=0)
SAMP = pd.read_table('1124_ALL_SAMP.list',header=None)
for SRR in SAMP[0]:
	MU1_W_COUNT[SRR] = 0
	W = pd.read_table(f'{SRR}_MU1.W_count.tsv', names=['name','W_COUNT'])
	for IL in W['name']:
		MU1_W_COUNT.loc[MU1_W_COUNT['name']==IL,SRR] = W.loc[W['name']==IL,'W_COUNT'].to_list()[0]
MASK = pd.read_table('MU1_chrAIL_SURR500_DVSB_MASKED.intersect.txt',header=None)[3].drop_duplicates()
ID_CONVER = pd.read_table(f'MU1_IL_ILSURROUNDING.tsv',names=['name','surr'])
ID_CONVER = ID_CONVER.loc[ID_CONVER['surr'].isin(MASK)]
for SRR in SAMP[0]:
	MU1_MASKW_COUNT[SRR] = 0
	MASKW = pd.read_table(f'{SRR}_MU1.MASK_W_count.tsv', names=['name','W_COUNT'])
	for IL in MASKW['name']:
		if IL in ID_CONVER['name']:
			MU1_MASKW_COUNT.loc[MU1_MASKW_COUNT['name']==IL,SRR] = MASKW.loc[MASKW['name']==IL,'W_COUNT'].to_list()[0]
MU1_W_COUNT = MU1_W_COUNT.set_index(['name','DISC2IL'])
MU1_MASKW_COUNT = MU1_MASKW_COUNT.set_index(['name','DISC2IL'])
MU1_WMERGE_COUNT = MU1_W_COUNT + MU1_MASKW_COUNT
MU1_WMERGE_COUNT['MASK'] = 0
MU1_WMERGE_COUNT = MU1_WMERGE_COUNT.reset_index()
MU1_WMERGE_COUNT.loc[MU1_WMERGE_COUNT['name'].isin(ID_CONVER['name']),'MASK'] = 1
MU1_WMERGE_COUNT.to_csv('MU1_WMERGE_COUNT.tsv',sep="\t",index=False,header=True)

MU2_W_COUNT = pd.read_table('MU2_GENOTYPE.tsv', header=0)
MU2_MASKW_COUNT = pd.read_table('MU2_GENOTYPE.tsv', header=0)
SAMP = pd.read_table('1124_ALL_SAMP.list',header=None)

for SRR in SAMP[0]:
	MU2_W_COUNT[SRR] = 0
	W = pd.read_table(f'{SRR}_MU2.W_count.tsv', names=['name','W_COUNT'])
	for IL in W['name']:
		MU2_W_COUNT.loc[MU2_W_COUNT['name']==IL,SRR] = W.loc[W['name']==IL,'W_COUNT'].to_list()[0]

MASK = pd.read_table('MU2_chrAIL_SURR500_DVSB_MASKED.intersect.txt',header=None)[3].drop_duplicates()
ID_CONVER = pd.read_table(f'MU2_IL_ILSURROUNDING.tsv',names=['name','surr'])
ID_CONVER = ID_CONVER.loc[ID_CONVER['surr'].isin(MASK)]

for SRR in SAMP[0]:
	MU2_MASKW_COUNT[SRR] = 0
	MASKW = pd.read_table(f'{SRR}_MU2.MASK_W_count.tsv', names=['name','W_COUNT'])
	for IL in MASKW['name']:
		if IL in ID_CONVER['name']:
			MU2_MASKW_COUNT.loc[MU2_MASKW_COUNT['name']==IL,SRR] = MASKW.loc[MASKW['name']==IL,'W_COUNT'].to_list()[0]

MU2_W_COUNT = MU2_W_COUNT.set_index(['name','DISC2IL'])
MU2_MASKW_COUNT = MU2_MASKW_COUNT.set_index(['name','DISC2IL'])
MU2_WMERGE_COUNT = MU2_W_COUNT + MU2_MASKW_COUNT
MU2_WMERGE_COUNT['MASK'] = 0
MU2_WMERGE_COUNT = MU2_WMERGE_COUNT.reset_index()
MU2_WMERGE_COUNT.loc[MU2_WMERGE_COUNT['name'].isin(ID_CONVER['name']),'MASK'] = 1
MU2_WMERGE_COUNT.to_csv('MU2_WMERGE_COUNT.tsv',sep="\t",index=False,header=True)
```

# Check all the ILs in all samples

```bash
for SAMP in $( cat 1124_ALL_SAMP.list ); do
echo ${SAMP} > MU1/${SAMP}.MU1.SPLIT_COUNT.txt
echo ${SAMP} > MU2/${SAMP}.MU2.SPLIT_COUNT.txt
echo ${SAMP} > MU1/${SAMP}.MU1.DISC_COUNT.txt
echo ${SAMP} > MU2/${SAMP}.MU2.DISC_COUNT.txt
	bedtools intersect -wa -c -a ALL_MU1.2SUP.bed -b MU1/${SAMP}.MU1.SPLIT_SUP.bed | \
cut -f 4 >> MU1/${SAMP}.MU1.SPLIT_COUNT.txt
	bedtools intersect -wa -c -a ALL_MU2.2SUP.bed -b MU2/${SAMP}.MU2.SPLIT_SUP.bed | \
cut -f 4 >> MU2/${SAMP}.MU2.SPLIT_COUNT.txt
	bedtools intersect -wa -c -a ALL_MU1.2SUP.bed -b MU1/${SAMP}.MU1.DISC_SUP.bed | \
cut -f 4 >> MU1/${SAMP}.MU1.DISC_COUNT.txt
	bedtools intersect -wa -c -a ALL_MU2.2SUP.bed -b MU2/${SAMP}.MU2.DISC_SUP.bed | \
cut -f 4 >> MU2/${SAMP}.MU2.DISC_COUNT.txt
done

cp ALL_MU1.DISC2IL_marked.tsv TEMP1
for SAMP in $( cat 1124_ALL_SAMP.list ); do
paste TEMP1 MU1/${SAMP}.MU1.SPLIT_COUNT.txt > TEMP2
mv TEMP2 TEMP1
done
mv TEMP1 ALL_MU1.DISC2IL_marked.244.SPLIT_COUNT.tsv

cp ALL_MU1.DISC2IL_marked.tsv TEMP1
for SAMP in $( cat 1124_ALL_SAMP.list ); do
paste TEMP1 MU1/${SAMP}.MU1.DISC_COUNT.txt > TEMP2
mv TEMP2 TEMP1
done
mv TEMP1 ALL_MU1.DISC2IL_marked.244.DISC_COUNT.tsv

cp ALL_MU2.DISC2IL_marked.tsv TEMP1
for SAMP in $( cat 1124_ALL_SAMP.list ); do
paste TEMP1 MU2/${SAMP}.MU2.SPLIT_COUNT.txt > TEMP2
mv TEMP2 TEMP1
done
mv TEMP1 ALL_MU2.DISC2IL_marked.244.SPLIT_COUNT.tsv

cp ALL_MU2.DISC2IL_marked.tsv TEMP1
for SAMP in $( cat 1124_ALL_SAMP.list ); do
paste TEMP1 MU2/${SAMP}.MU2.DISC_COUNT.txt > TEMP2
mv TEMP2 TEMP1
done
mv TEMP1 ALL_MU2.DISC2IL_marked.244.DISC_COUNT.tsv

```

## MULE activity in RNAseq data

### **#prepare the reference for RNAseq analysis**

```bash
cut -f 1,2,3 /zfs/socbd/bwu4/YU/annotation/FOUR_TE/TE_IN_CITRUS_GENOMES/CORRECTED_DVS_4TE.bed > CORRECTED_DVS_MULE.raw.bed

#>Manually add the 5 non-intact MULE regions in CORRECTED_DVS_4TE.1-3.bed

bedtools sort -i CORRECTED_DVS_MULE.raw.bed > CORRECTED_DVS_MULE.bed 

# Output all genic regions overlapping with any MULEs
bedtools intersect -wo -a CORRECTED_DVS_MULE.bed -b CK2021.FINAL0403.masked.gff3 | awk '$6=="gene" {print $4"\t"$7"\t"$8}' > MU_OVERLAP_GENE.bed
# Mask MULE overlapping genic regions in gff3 and the preivously masked genome for RNAseq analysis
bedtools intersect -v -a CK2021.FINAL0403.masked.gff3 -b MU_OVERLAP_GENE.bed > DVS.20211128.MASKED.gff3
bedtools maskfasta -fi CK2021.FINAL0403.masked.fasta -bed CORRECTED_DVS_MULE.bed -fo DVS.20211128.MASKED.fasta

# Get the modified representative MULE sequences
infoseq -only -name /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/CK2021.09282021.MASKED.fasta | grep Mu > Mu.id.txt
cdbfasta /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/CK2021.09282021.MASKED.fasta
cat Mu.id.txt | cdbyank /zfs/socbd/bwu4/YU/annotation/FOUR_TE/09212021_NEW/CK2021.09282021.MASKED.fasta.cidx - > Modified_Mu.fasta
bedtools maskfasta -fi CK2021.FINAL0403.masked.fasta -bed CORRECTED_DVS_MULE.bed -fo DVS.20211128.MASKED.fasta
cat Modified_Mu.fasta >> DVS.20211128.MASKED.fasta

# Index the reference for HISAT2
gffread -E DVS.20211128.MASKED.gff3 -T > DVS.20211128.MASKED.gtf
extract_splice_sites.py DVS.20211128.MASKED.gtf > DVS.20211128.MASKED.ss
extract_exons.py DVS.20211128.MASKED.gtf > DVS.20211128.MASKED.exon
hisat2-build --ss *.ss --exon *.exon DVS.20211128.MASKED.fasta DVS.20211128.MASKED

# move all files to /zfs/socbd/bwu4/MULE/RNAseq/REF
mv * /zfs/socbd/bwu4/MULE/RNAseq/REF
```

### #mapping RNAseq data to the modified reference genome

```bash
cd /scratch1/bwu4/SWO_RNASEQ/SAM
hisat2 -p 40 --dta -x /zfs/socbd/bwu4/MULE/RNAseq/REF/DVS.20211128.MASKED \
-1 ../SRR12096786_1.clean.fastq -2 ../SRR12096786_2.clean.fastq -S SRR12096786_DVS.20211128.MASKED.sam
samtools sort -@ 8 -o SRR12096786_DVS.20211128.MASKED.srt.bam SRR12096786_DVS.20211128.MASKED.sam 
samtools index -@ 8 SRR12096786_DVS.20211128.MASKED.srt.bam

hisat2 -p 40 --dta -x /zfs/socbd/bwu4/MULE/RNAseq/REF/DVS.20211128.MASKED \
-1 ../SRR12096786_1.clean.fastq -2 ../SRR12096786_2.clean.fastq -S SRR12096786_DVS.20211128.MASKED.sam
samtools sort -@ 8 -o SRR12096786_DVS.20211128.MASKED.srt.bam SRR12096786_DVS.20211128.MASKED.sam 
samtools index -@ 8 SRR12096786_DVS.20211128.MASKED.srt.bam

###### TEMPLATE.PBS ######
#PBS -N TEMP
#PBS -l select=1:ncpus=8:mem=30gb:interconnect=1g,walltime=24:00:00
#PBS -j oe
#PBS -o TEMP.HISAT2.log

source ~/.bashrc
cd /scratch1/bwu4/SWO_RNASEQ/SAM
hisat2 -p 8 --dta -x /zfs/socbd/bwu4/MULE/RNAseq/REF/DVS.20211128.MASKED \
-1 ../TEMP_1.clean.fastq -2 ../TEMP_2.clean.fastq -S TEMP_DVS.20211128.MASKED.sam

###### TEMPLATE.PBS ends here ######

mkdir /scratch1/bwu4/SWO_RNASEQ/PBS && cd /scratch1/bwu4/SWO_RNASEQ/PBS
for SRR in $( cat ../SWO_RNASEQ.SRR | head -n 740 ); do
sed "s/TEMP/${SRR}/g" TEMPLATE.PBS > ${SRR}.pbs
done

array=(SRR8186430 SRR12096786)
for SRR in $( cat ../SWO_RNASEQ.SRR | head -n 740 ); do
if [[ ! " ${array[*]} " =~ ${SRR} ]]; then 
N=$( qstat -u bwu4 | grep RR | wc -l)
while [ $N -ge 100 ]; do
sleep 60
N=$( qstat -u bwu4 | grep RR | wc -l)
done
qsub ${SRR}.pbs
fi
done

###### SORTBAM.PBS ######
#PBS -N TEMP
#PBS -l select=1:ncpus=2:mem=30gb:interconnect=1g,walltime=24:00:00
#PBS -j oe
#PBS -o TEMP.SORTBAM.log

source ~/.bashrc
cd /scratch1/bwu4/SWO_RNASEQ/BAM
samtools sort -o TEMP_DVS.20211128.MASKED.srt.bam ../SAM/TEMP_DVS.20211128.MASKED.sam 
samtools index TEMP_DVS.20211128.MASKED.srt.bam
###### SORTBAM.PBS ends here ######

for SRR in $( cat ../SWO_RNASEQ.SRR | head -n 740 ); do
sed "s/TEMP/${SRR}/g" SORTBAM.PBS > ../PBS/${SRR}.SORTBAM.pbs
done

for SRR in $( cat ../SWO_RNASEQ.SRR | head -n 740 ); do
if [[ ! ${SRR} == "SRR12096786" ]]; then 
N=$( qstat -u bwu4 | grep RR | wc -l)
while [ $N -ge 100 ]; do
sleep 60
N=$( qstat -u bwu4 | grep RR | wc -l)
done
qsub ${SRR}.SORTBAM.pbs
fi
done
```

### Calculate FPM of Mu1 and Mu2 for all the samples

```bash
# Output the IDs of read pairs overlapping with Mu1 and Mu2
for SAMP in $( cat 1124_ALL_SAMP.list ); do
samtools view -F 268 BAM/${SAMP}.TE_masked.srt.bam | awk '$3~/Mu1/ {print $1}' | sort -u > BAM/${SAMP}.Mu1.readid
samtools view -q 1 -F 268 BAM/${SAMP}.TE_masked.srt.bam | awk '$3~/Mu2/ {print $1}' | sort -u > BAM/${SAMP}.Mu2.readid
done

for SAMP in $( tail -n 241 1124_ALL_SAMP.list ); do samtools stat -@ 28 BAM/${SAMP}.TE_masked.srt.bam > BAM/${SAMP}.BAM_stat.report ; done

for SAMP in $( cat 1124_ALL_SAMP.list ); do
# The deleted region : Mu2_1:3234-4981 , only reads with primary alignment overlapping with Mu2_1:(3234+50)-(4981-50) were counted
samtools view -F 268 BAM/${SAMP}.TE_masked.srt.bam Mu2_1:3284-4931 | cut -f 1 | sort -u > BAM/${SAMP}.Mu2_DEL.readid
done

for SAMP in $( cat 1124_ALL_SAMP.list ); do
# COUNT is the number of read pairs (segments) overlapping with Mu2
COUNT=$( wc -l BAM/${SAMP}.Mu2.readid | cut -d " " -f 1)
# DEL is the number of read pairs (segments) overlapping with the deleted Mu2_1:3234-4981 region in DVS_Mu2_2
DEL=$( wc -l BAM/${SAMP}.Mu2_DEL.readid | cut -d " " -f 1 )
# MAPPED is the total number of mapped read pairs of ${SAMP} 
MAPPED=$( grep "reads mapped:" BAM/${SAMP}.BAM_stat.report | awk '{print $4/2}' )
echo -e "${SAMP}\t${COUNT}\t${DEL}\t${MAPPED}" >> DNA244_Mu2.readpair_count.tsv
done

for SAMP in $( cat 1124_ALL_SAMP.list ); do
# COUNT is the number of read pairs (segments) overlapping with Mu2
COUNT=$( wc -l BAM/${SAMP}.Mu1.readid | cut -d " " -f 1)
# MAPPED is the total number of mapped read pairs of ${SAMP} 
MAPPED=$( grep "reads mapped:" BAM/${SAMP}.BAM_stat.report | awk '{print $4/2}' )
echo -e "${SAMP}\t${COUNT}\t${MAPPED}" >> DNA244_Mu1.readpair_count.tsv
done

```

### **Calculate FPM of Mu2 in 740 sweet orange RNAseq data**

```bash
awk '$3~/rRNA/ {print $1"\t"$4-1"\t"$5}' /zfs/socbd/bwu4/YU/annotation/CK.03022021.rfam_infernal.gff3 > DVS.rRNA.bed

for SAMP in $( head -n 740 SWO_RNASEQ.SRR ); do
samtools view -F 256 BAM/${SAMP}_DVS.20211128.MASKED.srt.bam Mu2_1 | grep chr > ${SAMP}.Mu2_chr.SAM
COUNT=$( samtools view -F 256 BAM/${SAMP}_DVS.20211128.MASKED.srt.bam Mu2_1 | cut -f 1 | sort -u | wc -l )
TOTAL=$( awk '$0~/reads;/ {match($0,/^([0-9]+) reads/,a); print a[1]}' PBS/${SAMP}.HISAT2.log )
MAPPED=$( awk '$0~/overall/ {match($0,/^([0-9\.]+)%/,a); print a[1]}' PBS/${SAMP}.HISAT2.log )
echo -e "${SAMP}\t${COUNT}\t${TOTAL}\t${MAPPED}" >> RNA740_Mu2.count.tsv
done
for SAMP in $( head -n 740 SWO_RNASEQ.SRR ); do
rRNA_COUNT=$( samtools view -@ 20 -F 256 -L DVS.rRNA.bed BAM/${SAMP}_DVS.20211128.MASKED.srt.bam | cut -f 1 | sort -u | wc -l )
echo -e "${SAMP}\t${rRNA_COUNT}" >> RNA740_Mu2.rRNA_count.tsv
done

cd /scratch1/bwu4/MANDARIN_RNASEQ_12272021
for BAM in *.bam ; do echo ${BAM%_DVS.20211128.MASKED.srt.bam} >> ../MAN_RNASEQ39.SRR ; done
for SAMP in $( cat MAN_RNASEQ39.SRR ); do
samtools view -F 256 BAM/${SAMP}_DVS.20211128.MASKED.srt.bam Mu2_1 | grep chr > ${SAMP}.Mu2_chr.SAM
COUNT=$( samtools view -F 256 BAM/${SAMP}_DVS.20211128.MASKED.srt.bam Mu2_1 | cut -f 1 | sort -u | wc -l )
TOTAL=$( awk '$0~/reads;/ {match($0,/^([0-9]+) reads/,a); print a[1]}' PBS/${SAMP}.HISAT2.log )
MAPPED=$( awk '$0~/overall/ {match($0,/^([0-9\.]+)%/,a); print a[1]}' PBS/${SAMP}.HISAT2.log )
echo -e "${SAMP}\t${COUNT}\t${TOTAL}\t${MAPPED}" >> MAN_RNA39_Mu2.count.tsv
done
for SAMP in $( cat MAN_RNASEQ39.SRR ); do
rRNA_COUNT=$( samtools view -@ 20 -F 256 -L DVS.rRNA.bed BAM/${SAMP}_DVS.20211128.MASKED.srt.bam | cut -f 1 | sort -u | wc -l )
echo -e "${SAMP}\t${rRNA_COUNT}" >> MAN_RNA39.rRNA_count.tsv
done

cd /scratch1/bwu4/PUMMELO_RNASEQ_12272021/BAM
for BAM in *.bam ; do echo ${BAM%_DVS.20211128.MASKED.srt.bam} >> ../PUM_RNASEQ171.SRR ; done
cd ..
for SAMP in $( cat PUM_RNASEQ171.SRR ); do
samtools view -F 256 BAM/${SAMP}_DVS.20211128.MASKED.srt.bam Mu2_1 | grep chr > ${SAMP}.Mu2_chr.SAM
COUNT=$( samtools view -F 256 BAM/${SAMP}_DVS.20211128.MASKED.srt.bam Mu2_1 | cut -f 1 | sort -u | wc -l )
TOTAL=$( awk '$0~/reads;/ {match($0,/^([0-9]+) reads/,a); print a[1]}' PBS/${SAMP}.HISAT2.log )
MAPPED=$( awk '$0~/overall/ {match($0,/^([0-9\.]+)%/,a); print a[1]}' PBS/${SAMP}.HISAT2.log )
echo -e "${SAMP}\t${COUNT}\t${TOTAL}\t${MAPPED}" >> PUM_RNA171_Mu2.count.tsv
done
for SAMP in $( cat PUM_RNASEQ171.SRR ); do
rRNA_COUNT=$( samtools view -@ 20 -F 256 -L DVS.rRNA.bed BAM/${SAMP}_DVS.20211128.MASKED.srt.bam | cut -f 1 | sort -u | wc -l )
echo -e "${SAMP}\t${rRNA_COUNT}" >> PUM_RNA171.rRNA_count.tsv
done

```

### Get read count for analyzed genes

```bash

for SAMP in $( head -n 740 SWO_RNASEQ.SRR ); do
COUNT=()
for i in {1..60}; do
input=($( sed -n ${i}p gene.list ))
COUNT=(echo -e ${COUNT}"\t"$( samtools view -q 1 BAM/${SAMP}_DVS.20211128.MASKED.srt.bam ${input[1]} | cut -f 1 | sort -u | wc -l ))
done
echo -e "${SAMP}\t${COUNT}" >> RNA740_GENES.count.tsv

for SAMP in $( head -n 740 SWO_RNASEQ.SRR ); do
rRNA_COUNT=$( samtools view -@ 20 -F 256 -L DVS.rRNA.bed BAM/${SAMP}_DVS.20211128.MASKED.srt.bam | cut -f 1 | sort -u | wc -l )
echo -e "${SAMP}\t${rRNA_COUNT}" >> RNA740_Mu2.rRNA_count.tsv
done

for SAMP in $( head -n 740 SWO_RNASEQ.SRR ); do 
 COUNT=() 
for i in {1..60}; do
 input=($( sed -n ${i}p gene.list )) 
 COUNT=$(echo -e ${COUNT}"\t"$( samtools view -@ 10 -q 1 BAM/${SAMP}_DVS.20211128.MASKED.srt.bam ${input[1]} | cut -f 1 | sort -u | wc -l )) 
done 
 echo -e "${SAMP}\t${COUNT}" >> RNA740_GENES.count.tsv 
done

for SAMP in $( head -n 740 SWO_RNASEQ.SRR ); do 
 COUNT=$( samtools view -@ 10 -q 1 BAM/${SAMP}_DVS.20211128.MASKED.srt.bam chr6B:24509960-24518502 | cut -f 1 | sort -u | wc -l ) 
 echo -e "${SAMP}\t${COUNT}" >> RNA740_DSWO6B02241.count.tsv 
done

```

### Output 3’ downstream terminal of transcribed MU2 mem

```bash
MU2_IL.bed
for SAMP in $( head -n 740 SWO_RNASEQ.SRR ); do
awk '$4>6500 {print $7"\t"$8-1"\t"$8}' SWO740_MU2_CHR/${SAMP}.Mu2_chr.SAM  | bedtools sort -i - | bedtools intersect -c -a MU2_IL.bed -b - > MU2_CHR_${SAMP}.COUNT.bed
done

cut -f 1,2,3,4 MU2_CHR_SRR8186430.COUNT.bed > TEMP
for SAMP in $( head -n 740 SWO_RNASEQ.SRR ); do
cut -f 5 MU2_CHR_${SAMP}.COUNT.bed | paste TEMP - > TEMP1
mv TEMP1 TEMP
done
mv TEMP SWO740_MU2_CHR.COUNT.bed
```
