# Detection of active transposable elements (TE) and TE insertion locus (IL) scanning

## Software requirements
 - Python3.* (version>=3.7; tested on 3.9 and 3.10)
 - gawk (tested on GNU Awk 5.1.0)
 - cd-hit (version=4.8.1; only used for active TE family detection)
 - samtools (version>=1.12; tested on v1.15, v1.17, and v1.18)
 - bcftools (version>=1.12; tested on v1.15, v1.17, and v1.18)
 - blast (version>=2.14; tested on v2.14)
 - pandas (tested on v1.5.2)
 - bwa (tested on version=0.7.17)
 - emboss (tested on version=6.6.0)
 - minimap2 (version>=2.17; tested on v2.17 and v2.26)
    
## One-step software installation
```bash
# Download the repository
git clone https://github.com/TheLuoFengLab/Transposon-insertion-locus-simulation-and-scanning.git

# Install tools in TESCAN.yml
conda env create -n TESCAN -f Transposon-insertion-locus-simulation-and-scanning/Conda_environment/TESCAN.yml

# Or, install using TESCAN.txt
conda create --name TESCAN --file Transposon-insertion-locus-simulation-and-scanning/Conda_environment/TESCAN.txt

# Activate conda environment
conda activate TESCAN
```
## 1. Active TE in SWO assemblies
See 
1.1 Mapping 10 SWO assemblies to DVS and call variants

```bash
for FAS in GCA_019144245.1.fasta GCA_019144225.1.fasta GCA_019144195.1.fasta GCA_019144185.1.fasta GCA_019144155.1.fasta GCA_019143665.1.fasta GCA_018104345.1.fasta Csiv4.fasta T78.asem.fasta SF.asem.fasta ; do
minimap2 -cx asm5 -t8 --cs CK2021.60.corrected.fasta $FAS > asm.paf  
sort -k6,6 -k8,8n asm.paf > asm.srt.paf             # sort by reference start coordinate
k8 paftools.js call asm.srt.paf > ${FAS}.var.txt
done

for VAR in *.fasta.var.txt ; do
awk '$5==1 && $4-$3>50 && $4-$3<20000 && $11-$10<=10 {print $2"\t"$3"\t"$4}' $VAR > ${VAR%.fasta.var.txt}.DEL.tsv
awk '$5==1 && $4-$3<=10 && $11-$10>50 && $11-$10<20000 {print $2"\t"$3"\t"$4}' $VAR > ${VAR%.fasta.var.txt}.INS.tsv

cat ${VAR%.fasta.var.txt}.DEL.tsv ${VAR%.fasta.var.txt}.INS.tsv > ${VAR%.fasta.var.txt}.INDEL.bed
done

for BED in Csiv4.INDEL.bed T78.asem.INDEL.bed SF.asem.INDEL.bed GCA_019144245.1.INDEL.bed GCA_019144225.1.INDEL.bed GCA_019144195.1.INDEL.bed GCA_019144185.1.INDEL.bed GCA_019144155.1.INDEL.bed GCA_019143665.1.INDEL.bed GCA_018104345.1.INDEL.bed; do
echo $BED
bedtools intersect -f 0.95 -F 0.95 -c -a ALL10.INDEL.bed -b $BED > ALL10.${BED}
done

for BED in ALL10.*.INDEL.bed ; do cut -f 4 $BED | paste TEMP - > TEMP1 ; mv TEMP1 TEMP; done

paste ALL10.INDEL.bed TEMP | sed "s/\t\t/\t/g" > TEMP1
mv TEMP1 TEMP
echo -e "CHR\tSTART\tEND\t"$( ls ALL10.*.INDEL.bed ) | sed "s/ /\t/g" > header.txt

cat header.txt TEMP > ALL10.INDEL.tsv
# Remove duplicate indels with both ends with distances <= 15 bp in excel
# Name all INDELs in excel and paste the deletions in ALL10.INS.bed
sed -i 's/INDEL/DEL/g' ALL10.DEL.bed
bedtools getfasta -fi CK2021.60.corrected.fasta -bed ALL10.DEL.bed -fo ALL10.DEL.fasta

awk '$4-$3<=10 && $5==1 && $11-$10>=50' *.var.txt | sort -u -k1,3 | awk '{print ">"$2":"$3"-"$4"\t"$8}' | sed "s/\t/\n/g" > INS.fasta

cat DEL.gi | cdbyank INS.fasta.cidx > ALL10.INS.fasta

cat ALL10.INS.fasta ALL10.DEL.fasta > ALL10.INDEL.fasta

cd-hit-est -r 1 -g 1 -c 0.80 -i ALL10.INDEL.fasta -o ALL10.INDEL.80_50.clusters.fasta -T 0 -aL 0.50 -M 370000 -d 50 -n 10
awk '$0~/>Cluster/ {print a"\t"b"\t"c"\t"e; a=$0;b=0;c=0} $0!~/>Cluster/ {b+=1; if ($0~/*/) {match($0,/([0-9]+)nt, >(.+)\.\.\./,d);c=d[2];e=d[1]}} END {print a"\t"b"\t"c"\t"e}' ALL10.INDEL.80_50.clusters.fasta.clstr | sed "s/ //g" > ALL10.INDEL.80_50.clusters.stat.tsv

cd-hit-est -r 1 -g 1 -c 0.95 -i ALL10.INDEL.fasta -o ALL10.INDEL.95_90.clusters.fasta -T 0 -aL 0.90 -M 370000 -d 50 -n 10
awk '$0~/>Cluster/ {print a"\t"b"\t"c"\t"e; a=$0;b=0;c=0} $0!~/>Cluster/ {b+=1; if ($0~/*/) {match($0,/([0-9]+)nt, >(.+)\.\.\./,d);c=d[2];e=d[1]}} END {print a"\t"b"\t"c"\t"e}' ALL10.INDEL.95_90.clusters.fasta.clstr | sed "s/ //g" > ALL10.INDEL.95_90.clusters.stat.tsv

```


#Maually curate all TE member terminal haplotypes from 32 different families and output all of them as NEW*_LEFT/RIGHT_50bp.align.fasta

#Python script for summarizing terminal haplotypes d:\SWO\SWO_TE\TERMINAL_HAP.py

```python
#!/usr/bin/env python3

import pandas as pd
import os
import sys

PRE=sys.argv[1]

RIGHT_ALL=pd.read_table(PRE+'_RIGHT_50bp.align.fasta',names=['ID','R_HAP','R_COUNT']).drop(['R_COUNT'],axis=1)
RIGHT_COUNT = RIGHT_ALL.groupby('R_HAP').count().reset_index().rename(columns={"ID": "R_COUNT"})
RIGHT_COUNT = RIGHT_COUNT.sort_values('R_COUNT',ascending=False)
RIGHT_COUNT = RIGHT_COUNT.reset_index(drop=True)
RIGHT_COUNT['RHAP_ID'] = 'R'+RIGHT_COUNT.index.astype(str)
RIGHT_ALL = RIGHT_ALL.merge(RIGHT_COUNT, on='R_HAP')

LEFT_ALL=pd.read_table(PRE+'_LEFT_50bp.align.fasta',names=['ID','L_HAP','L_COUNT']).drop(['L_COUNT'],axis=1)
LEFT_COUNT = LEFT_ALL.groupby('L_HAP').count().reset_index().rename(columns={"ID": "L_COUNT"})
LEFT_COUNT = LEFT_COUNT.sort_values('L_COUNT',ascending=False)
LEFT_COUNT = LEFT_COUNT.reset_index(drop=True)
LEFT_COUNT['LHAP_ID'] = 'L'+LEFT_COUNT.index.astype(str)
LEFT_ALL = LEFT_ALL.merge(LEFT_COUNT, on='L_HAP')

ALL = LEFT_ALL.merge(RIGHT_ALL,on='ID')
ALL.to_csv(PRE+'_MEMBER_TERMINAL_HAP.tsv', sep="\t", header=True, index=False)
ALL_SUB=ALL.loc[(ALL['L_COUNT']>1) | (ALL['R_COUNT']>1) ]

LEFT=ALL_SUB[['LHAP_ID', 'L_HAP']].drop_duplicates()
LEFT['LHAP_ID']=">"+LEFT['LHAP_ID']
LEFT.to_csv(PRE+'_LEFT_50bp_SUB.fasta', sep="\t", header=False, index=False)

RIGHT=ALL_SUB[['RHAP_ID', 'R_HAP']].drop_duplicates()
RIGHT['RHAP_ID']=">"+RIGHT['RHAP_ID']
RIGHT.to_csv(PRE+'_RIGHT_50bp_SUB.fasta', sep="\t", header=False, index=False)
```

```bash
cd /mnt/d/SWO/SWO_TE
for N in {1..32} ; do x="NEW"$N ; ./TERMINAL_HAP.py $x ; done
```

```bash
for N in {1..32} ; do 
blastn -task blastn-short -num_threads 40 -query NEW${N}_LEFT_50bp_SUB.fasta \
-db CK2021 -outfmt 6 | awk '$3==100 && $4==50 && $9<$10 {print $1"\t"$2"\t"$9"\tPLUS"} $3==100 && $4==50 && $9>$10 {print $1"\t"$2"\t"$9"\tMINUS"}' > NEW${N}_PRE_LEFT.tsv
blastn -task blastn-short -num_threads 40 -query NEW${N}_RIGHT_50bp_SUB.fasta \
-db CK2021 -outfmt 6 | awk '$3==100 && $4==50 && $9<$10 {print $1"\t"$2"\t"$10"\tPLUS"} $3==100 && $4==50 && $9>$10 {print $1"\t"$2"\t"$10"\tMINUS"}' > NEW${N}_PRE_RIGHT.tsv
done
```

```bash
for ID in {2..32}; do ./IDENTIFY_COPY_IN_DVS.py NEW${ID} ; done
```

### ### IDENTIFY_COPY_IN_DVS.py

```python
#!/usr/bin/env python3

import pandas as pd
import os
import sys

PRE=sys.argv[1]

LEFT = pd.read_table(PRE+'_PRE_LEFT.tsv',names=['LHAP_ID','CHR','L_COOR','STRAND'])
RIGHT = pd.read_table(PRE+'_PRE_RIGHT.tsv',names=['RHAP_ID','CHR','R_COOR','STRAND'])

LIST=[]
for i in LEFT.index:
    if LEFT.iloc[i]['STRAND']=="MINUS" :
        DF=RIGHT.loc[(RIGHT['CHR']==LEFT.iloc[i]['CHR']) & (LEFT.iloc[i]['L_COOR']-RIGHT['R_COOR']<20000) & (LEFT.iloc[i]['L_COOR']-RIGHT['R_COOR']>0) & (LEFT.iloc[i]['STRAND']==RIGHT['STRAND'])]
    if LEFT.iloc[i]['STRAND']=="PLUS" :
        DF=RIGHT.loc[(RIGHT['CHR']==LEFT.iloc[i]['CHR']) & (LEFT.iloc[i]['L_COOR']-RIGHT['R_COOR']>-20000) & (LEFT.iloc[i]['L_COOR']-RIGHT['R_COOR']<0) & (LEFT.iloc[i]['STRAND']==RIGHT['STRAND'])]
    if len(DF.index)>0:
        DF['LHAP_ID']=LEFT.iloc[i]['LHAP_ID']
        DF['L_COOR']=LEFT.iloc[i]['L_COOR']
    else:
        DF=pd.DataFrame(LEFT.iloc[i]).T
        DF['RHAP_ID']="NA"
        DF['R_COOR']="NA"
    LIST.append(DF)

df=pd.concat(LIST)

df1=RIGHT.loc[~RIGHT['R_COOR'].isin(df['R_COOR'])]
df1['LHAP_ID']="NA"
df1['L_COOR']="NA"

pd.concat([df,df1]).to_csv(PRE+'_COPY_IN_DVS.tsv',sep="\t",header=True,index=False)
```

```bash
for N in {1..32} ; do
cat NEW${N}_COPY_IN_DVS.tsv | awk \
'NR>1 && $3!="NA" && $6!="NA" && $3<$6 {print $2"\t"$3-1"\t"$6"\t"$2":"$3-1"-"$6"\t.\t-"} NR>1 && $3!="NA" && $6!="NA" && $3>$6 {print $2"\t"$6-1"\t"$3"\t"$2":"$6-1"-"$3"\t.\t+"}' > NEW${N}_INTACT_COPIES.bed
cat NEW${N}_COPY_IN_DVS.tsv | awk \
'$3=="NA" && $4=="PLUS" {print $2"\t"$6-1"\t"$6+49} $3=="NA" && $4=="MINUS" {print $2"\t"$6-50"\t"$6} $6=="NA" && $4=="MINUS" {print $2"\t"$3-1"\t"$3+49} $6=="NA" && $4=="PLUS" {print $2"\t"$3-50"\t"$3}' > NEW${N}_UNINTACT_COPIES.bed
done

```

```bash
bedtools getfasta -s -fi CK2021.60.corrected.fasta -bed NEW8_INTACT_COPIES.bed  -fo NEW8_INTACT_COPIES.fasta
cd-hit-est -r 1 -g 1 -c 0.97 -i NEW8_INTACT_COPIES.fasta -o NEW8_INTACT_COPIES.97_95.clusters.fasta -T 0 -aL 0.95 -M 370000 -d 50 -n 10
```

### # REP_INTACT_BLAST.sh

```bash
#!/bin/bash
ID=$1
bedtools getfasta -s -fi CK2021.60.corrected.fasta -bed ${ID}_INTACT_COPIES.bed  -fo ${ID}_INTACT_COPIES.fasta
sed -i 's/(+)//g;s/(-)//g' ${ID}_INTACT_COPIES.fasta
cdbfasta ${ID}_INTACT_COPIES.fasta
grep ">" ${ID}_INTACT_COPIES.fasta | sort -u | sed "s/>//g" > ${ID}_INTACT_UNIQ.gi
cat ${ID}_INTACT_UNIQ.gi | cdbyank ${ID}_INTACT_COPIES.fasta.cidx - > ${ID}_INTACT_UNIQ_COPIES.fasta
makeblastdb -in ${ID}_INTACT_UNIQ_COPIES.fasta -out ${ID}_INTACT -parse_seqids -dbtype nucl
cd-hit-est -r 1 -g 1 -c 0.97 -i ${ID}_INTACT_COPIES.fasta -o ${ID}_INTACT_COPIES.97_95.clusters.fasta -T 0 -M 370000 -d 50 -n 10

cat ${ID}_REP.gi | cdbyank ${ID}_INTACT_COPIES.fasta.cidx - > ${ID}_REP.fasta
makeblastdb -in ${ID}_REP.fasta -out ${ID}_REP -parse_seqids -dbtype nucl
blastn -query ${ID}_REP.fasta -db ${ID}_REP -outfmt 6 -num_threads 40
blastn -query ${ID}_REP.fasta -db ${ID}_INTACT -outfmt 6 -num_threads 40 | awk -v x=$2 '$3>97 && $4>x' > ${ID}_REP_INTACT.results
```

### # TE representatives were selected for each TE family:

(1) Each terminal haplotype only selected once;

(2) Each cluster (${ID}_INTACT_COPIES.97_95.clusters.fasta.clstr) with ≥4 members has one member selected;

```bash
# Further check the uniqueness of the terminals of the TE representatives
for FASTA in *_50bp_SUB.fasta; do
blastn -task blastn-short -num_threads 40 -query ALLTE_TERMINAL.fasta -db ALLTE_REP -outfmt 6 | \
awk '$3==100 && $4==50' > ${FASTA%_50bp_SUB.fasta}_REP.results
done

cat NEW*_REP.fasta > ALLTE_REP.fasta
blastn -num_threads 40 -query ALLTE_REP.fasta -db CK2021 -outfmt 6 | \
awk '$3>95 && $4>=100 && $9<$10 {print $2"\t"$9-1"\t"$10} $3>97 && $4>=100 && $9>$10 {print $2"\t"$10-1"\t"$9}' > TO_MASK.bed
for TSV in NEW*_PRE_LEFT.tsv ; do
awk '$4=="PLUS" {print $2"\t"$3-1"\t"$3+49} $4=="MINUS" {print $2"\t"$3-50"\t"$3}' $TSV >> TE_LEFT.bed
done
for TSV in NEW*_PRE_RIGHT.tsv ; do
awk '$4=="PLUS" {print $2"\t"$3-50"\t"$3} $4=="MINUS" {print $2"\t"$3-1"\t"$3+49}' $TSV >> TE_RIGHT.bed
done

cat TO_MASK.bed TE_LEFT.bed TE_RIGHT.bed | bedtools sort -i - | bedtools merge -i - > NEW_TO_MASK.bed

blastn -query ALLTE_REP.fasta -db CK2021.05072022.MASKED -outfmt 6 -num_threads 40 | \
awk '$3>99 && $4>48 && $4<100 && $7<10 && $9<$10 {print $2"\t"$9-1"\t"$10} $3>99 && $4>48 && $4<100 && $7<10 && $9>$10 {print $2"\t"$10-1"\t"$9}' >> more_mask.bed

```

```bash
for SRR in $( head -n 10 SWO.list ); do
if [[ ! -a ${SRR}.TE_05082022.srt.bam ]]; then
bwa mem -t 40 CK2021.05072022.MASKED.fasta ${SRR}_1.clean.fastq \
${SRR}_2.clean.fastq > ${SRR}.sam
samtools sort -@ 20 -O BAM -o ${SRR}.TE_05082022.srt.bam ${SRR}.sam
fi
done
```

```bash
for SRR in $( cat LEN150.txt ); do
    if [[ ! -a ${SRR}.DISC_SUP.bed ]]; then  
        ./TE_INSERT_READ.sh TE_LIST.txt ${SRR} 40 150 500 TE CK2021.05072022.MASKED.fasta 
    fi 
done
```

## IL results statistics

```bash
# Separate MU1 and MU2 ILs
for pat in TE{1..32}; do
mkdir $pat
done

for SAMP in $( cat ALL0510.txt ); do
for pat in TE{1..32}_; do
awk -v pat=$pat '$4~pat && $2>0 {print $0}' ${SAMP}.SPLIT_SUP.bed > ${pat%_}/${SAMP}.${pat%_}.SPLIT_SUP.bed
done
done

for SAMP in $( cat ALL0510.txt ); do
for pat in TE{1..32}_; do
awk -v pat=$pat '$4~pat && $2>0 {print $0}' ${SAMP}.DISC_SUP.bed > ${pat%_}/${SAMP}.${pat%_}.DISC_SUP.bed
done
done

# Put all ILs detected by SPLIT and DISC reads in one bed file
for pat in TE{1..32}; do
cat ${pat}/*.${pat}.SPLIT_SUP.bed > ALL_${pat}.SPLIT.bed
done

for pat in TE{1..32}; do
cat ${pat}/*.${pat}.DISC_SUP.bed > ALL_${pat}.DISC.bed
done

for pat in TE{1..32}; do
    cat ALL_${pat}.SPLIT.bed ALL_${pat}.DISC.bed > ALL_${pat}.SPLIT_DISC.bed
    bedtools sort -i ALL_${pat}.SPLIT.bed | \
        bedtools merge -d 15 -i - > ALL_${pat}.SPLIT_IL.merged.bed
    bedtools intersect -wo -a ALL_${pat}.SPLIT_IL.merged.bed  \
        -b ALL_${pat}.SPLIT.bed > ALL_${pat}.SPLIT_IL.SPLIT.intersect
    bedtools intersect -wo -a ALL_${pat}.SPLIT_IL.merged.bed  \
        -b ALL_${pat}.DISC.bed > ALL_${pat}.SPLIT_IL.DISC.intersect
done

for pat in TE{1..32}; do
    bedtools intersect -c -a ALL_${pat}.SPLIT_IL.merged.bed \
        -b ALL_${pat}.SPLIT_DISC.bed > ALL_${pat}.SPLIT_IL.SPLIT_DISC.count.tsv
    awk '$4>1' ALL_${pat}.SPLIT_IL.SPLIT_DISC.count.tsv | \
        cut -f 1,2,3 > ALL_${pat}.2SUP.bed
done

for pat in TE{1..32}; do
    bedtools intersect -c -a ALL_${pat}.SPLIT_IL.merged.bed \
        -b ALL_${pat}.SPLIT.bed > ALL_${pat}.SPLIT_IL.SPLIT.count.tsv
    awk '$4>1' ALL_${pat}.SPLIT_IL.SPLIT.count.tsv | \
        cut -f 1,2,3 > ALL_${pat}.2SPLIT_SUP.bed
done

tr '\n' '\t' < ALL0510.txt > header.txt
echo -e "chr\tLeft\tRight\t$( cat header.txt )" > header.txt

bedtools intersect -c -a ALL_TE1.2SPLIT_SUP.bed -b TE1/CK.TE1.SPLIT_SUP.bed

for pat in TE{1..32}; do
rm -f TEMP.tsv
cp ALL_${pat}.2SPLIT_SUP.bed TEMP.tsv
for SAMP in $( cat ALL0510.txt ); do
bedtools intersect -c -a ALL_${pat}.2SPLIT_SUP.bed -b ${pat}/${SAMP}.${pat}.SPLIT_SUP.bed | \
cut -f 4 | paste TEMP.tsv - > TEMP1.tsv
mv TEMP1.tsv TEMP.tsv
done
cat header.txt TEMP.tsv > ALL_SAMP.${pat}.SPLIT.read_count.tsv
done

for pat in TE{1..32}; do
for SAMP in $( cat ALL0510.txt ); do
echo ${pat}" "${SAMP}"\t"$( wc -l ${pat}/${SAMP}.${pat}.IL.DISC_COUNT.bed )
done
done
```

 ## DISC_2IL.sh: Deal with DISC reads supporting more than one IL

```bash
#!/bin/bash 
pat=$1
for SAMP in $( cat ALL0510.txt ); do
# Identify ILs supported by split-mapped reads 
bedtools intersect -wa -a ALL_${pat}.2SPLIT_SUP.bed \
-b ${pat}/${SAMP}.${pat}.SPLIT_SUP.bed | \
sort -u > ${pat}/${SAMP}.${pat}_IL_SPLIT.bed
  
# Identify ILs not supported by split-mapped reads
bedtools intersect -v -a ALL_${pat}.2SPLIT_SUP.bed \
-b ${pat}/${SAMP}.${pat}.SPLIT_SUP.bed > ${pat}/${SAMP}.${pat}_IL_NOSPLIT.bed
  
# DISC reads uniquely supporting split-read supported ILs
bedtools intersect -wo -a ${pat}/${SAMP}.${pat}_IL_SPLIT.bed \
-b ${pat}/${SAMP}.${pat}.DISC_SUP.bed | cut -f 10 | sort | \
uniq -u > ${pat}/${SAMP}.${pat}_IL_SPLIT.uniq_DISC.txt

# ALL DISC (including those matching two ILs) reads supporting split-read supported ILs
bedtools intersect -wo -a ${pat}/${SAMP}.${pat}_IL_SPLIT.bed \
-b ${pat}/${SAMP}.${pat}.DISC_SUP.bed | cut -f 10 | \
sort -u > ${pat}/${SAMP}.${pat}_IL_SPLIT.ALL_DISC.txt

./SUBSET_COLUMN1.py ${pat} ${SAMP}

# DISC read count on split-read supported IL in SAMP
bedtools intersect -c -a ${pat}/${SAMP}.${pat}_IL_SPLIT.bed \
-b ${pat}/${SAMP}.${pat}.SPLIT_UNIQ_DISC_SUP.bed > \
${pat}/${SAMP}.${pat}.SPLIT_IL.DISC_COUNT.bed

# Output DISC read ids on ILs without split-read support
bedtools intersect -wo -a ${pat}/${SAMP}.${pat}_IL_NOSPLIT.bed \
-b ${pat}/${SAMP}.${pat}.NOSPLIT_DISC_SUP.bed | \
cut -f 10 | sort | uniq -u > ${SAMP}.${pat}.UNIQ_NOSPLIT_DISC_SUP.txt

./SUBSET_COLUMN2.py ${pat} ${SAMP}

# DISC read count on non-split-read supported IL in SAMP
bedtools intersect -c -a ${pat}/${SAMP}.${pat}_IL_NOSPLIT.bed \
-b ${pat}/${SAMP}.${pat}.NOSPLIT_UNIQ_DISC_SUP.bed > \
${pat}/${SAMP}.${pat}.NOSPLIT_IL.DISC_COUNT.bed

# Concatenate DISC read counts of SPLIT and NOSPLIT IL
cat ${pat}/${SAMP}.${pat}.SPLIT_IL.DISC_COUNT.bed ${pat}/${SAMP}.${pat}.NOSPLIT_IL.DISC_COUNT.bed | \
bedtools sort -i - > ${pat}/${SAMP}.${pat}.IL.DISC_COUNT.bed
done
```

```bash
parallel ./DISC_2IL.sh ::: TE{1..32}
```

### # SUBSET_COLUMN1.py

```python
#!/usr/bin/env python3

import pandas as pd
import sys

TE=sys.argv[1]
SAMP=sys.argv[2]

with open(f'{TE}/{SAMP}.{TE}_IL_SPLIT.uniq_DISC.txt') as file:
    UNIQ = file.readlines()
    UNIQ = [line.rstrip() for line in UNIQ]
with open(f'{TE}/{SAMP}.{TE}_IL_SPLIT.ALL_DISC.txt') as file:
    ALL = file.readlines()
    ALL = [line.rstrip() for line in ALL]

DISC=pd.read_table(f'{TE}/{SAMP}.{TE}.DISC_SUP.bed',header=None)
DISC.loc[DISC[6].isin(UNIQ)].to_csv(f'{TE}/{SAMP}.{TE}.SPLIT_UNIQ_DISC_SUP.bed', sep="\t", header=False, index=False)
DISC.loc[~DISC[6].isin(ALL)].to_csv(f'{TE}/{SAMP}.{TE}.NOSPLIT_DISC_SUP.bed', sep="\t", header=False, index=False)
```

## # SUBSET_COLUMN2.py

```python
#!/usr/bin/env python3

import pandas as pd
import sys

TE=sys.argv[1]
SAMP=sys.argv[2]

with open(f'{SAMP}.{TE}.UNIQ_NOSPLIT_DISC_SUP.txt') as file:
    UNIQ = file.readlines()
    UNIQ = [line.rstrip() for line in UNIQ]

DISC=pd.read_table(f'{TE}/{SAMP}.{TE}.DISC_SUP.bed', header=None)
DISC.loc[DISC[6].isin(UNIQ)].to_csv(f'{TE}/{SAMP}.{TE}.NOSPLIT_UNIQ_DISC_SUP.bed', sep="\t", header=False, index=False)
```

```bash
# Put all sample DISC read counts into one large table
for pat in TE{1..32}; do
rm -f TEMP.tsv
bedtools sort -i ALL_${pat}.2SPLIT_SUP.bed > TEMP.tsv
for SAMP in $( cat ALL0510.txt ); do
cut -f 4 ${pat}/${SAMP}.${pat}.IL.DISC_COUNT.bed | paste TEMP.tsv - > TEMP1.tsv
mv TEMP1.tsv TEMP.tsv
done
cat header.txt TEMP.tsv > ALL_SAMP.${pat}.DISC.read_count.tsv
done
```

 # Correspond NGS ILs with ILs in DVS 

```bash
for N in {1..32}; do

# Intact DVS TE member : chr1A:start-end
bedtools sort -i NEW${N}_INTACT_COPIES.EXT10.bed | bedtools merge -i - | \
awk '{print $1"\t"$2"\t"$3"\t"$1":"$2+10"-"$3-10}' > TE${N}_INTACT.EXT10.bed

# Unintact DVS TE member : chr1A|start-end
bedtools intersect -v -a NEW${N}_UNINTACT_COPIES.bed -b TE${N}_INTACT.EXT10.bed **| \**
awk '{print $1"\t"$2-10"\t"$3+10"\t"$1"|"$2"-"$3}' > TE${N}_UNINTACT_EXT10.bed

cat TE${N}_INTACT.EXT10.bed TE${N}_UNINTACT_EXT10.bed | \
bedtools sort -i - > TE${N}_DVS.INTACT_UNINTACT.bed 

bedtools intersect -wao -a ALL_TE${N}.2SPLIT_SUP.bed -b TE${N}_DVS.INTACT_UNINTACT.bed | \
awk '$7=="." {print $1"\t"$2"\t"$3"\tNOT_DVS"} $7!="." {print $1"\t"$2"\t"$3"\t"$7}' \
> ALL_TE${N}.2SPLIT_SUP.DVS.bed

done

for TE in TE{1..32}; do
awk '{print $1":"$2"-"$3"\t"$4}' ALL_${TE}.2SPLIT_SUP.DVS.bed > ALL_${TE}.2SPLIT_SUP.DVS.tsv
done

```

### # Content of SUM_SPLIT_DISC_COUNT.py

```python
#!/usr/bin/env python3

import pandas as pd
import sys

TE=sys.argv[1]
SPLIT=pd.read_table(f'ALL_SAMP.{TE}.SPLIT.read_count.tsv', header=0)
SPLIT.drop('Unnamed: 232', axis=1, inplace=True)
DISC=pd.read_table(f'ALL_SAMP.{TE}.DISC.read_count.tsv', header=0)
DISC.drop('Unnamed: 232', axis=1, inplace=True)

SPLIT.set_index(['chr', 'Left', 'Right'], inplace=True)
DISC.set_index(['chr','Left','Right'], inplace=True)

SUM=SPLIT+DISC
SUM.reset_index().to_csv(f'ALL_SAMP.{TE}.SPLIT_DISC_SUM.read_count.tsv', sep="\t", header=True, index=False)

SUM.where(SUM > 1, 0, inplace=True)
SUM.where(SUM < 2, 1, inplace=True)
SPLIT.where(SPLIT <1, 1, inplace=True)
GENOTYPE = SUM + SPLIT
GENOTYPE.where(GENOTYPE > 1, 0, inplace=True)
GENOTYPE.where(GENOTYPE < 2, 1, inplace=True)
GENOTYPE = GENOTYPE.reset_index()
GENOTYPE['name']=GENOTYPE['chr']+":"+GENOTYPE['Left'].astype(str)+"-"+GENOTYPE['Right'].astype(str)
GENOTYPE.drop(['chr', 'Left', 'Right'], axis = 1, inplace = True)
GENOTYPE.to_csv(f'{TE}_NEW_GENOTYPE.tsv', sep="\t", header=True, index=False)
```

### # Sum split and disc read count matrix

```bash
for TE in TE{1..32} ; do
./SUM_SPLIT_DISC_COUNT.py $TE
done
```

### Mutant-type cell ratio analysis

```bash
# Output all TE ILs in the DVS reference genome
for N in {1..32} ; do
awk '$4=="PLUS" {print $2"\t"$3-11"\t"$3-1} $4=="MINUS" {print $2"\t"$3"\t"$3+10}' /zfs/socbd/bwu4/YU/annotation/FOUR_TE/TE_IN_CITRUS_GENOMES/SWO/NEW${N}_PRE_LEFT.tsv >> DVS_TE${N}.bed
awk '$4=="PLUS" {print $2"\t"$3"\t"$3+10} $4=="MINUS" {print $2"\t"$3-11"\t"$3-1}' /zfs/socbd/bwu4/YU/annotation/FOUR_TE/TE_IN_CITRUS_GENOMES/SWO/NEW${N}_PRE_RIGHT.tsv >> DVS_TE${N}.bed
done
```

 # Calculate wild-type reads mapped on ILs

```bash
for TE in TE{1..32} ; do
for SAMP in $( cat ALL0510.txt ); do
rm -f ${SAMP}_${TE}.W_count.tsv
for IL in $( cat ${SAMP}_${TE}_IL.txt ); do
ILEXT=$( echo $IL | awk '{match($0,/(chr[1-9AB]+):([0-9]+)-([0-9]+)$/,a); print a[1]":"a[2]-1001"-"a[3]+1000}' )
W=$( samtools view -q 1 -f 66 ${SAMP}.TE_05082022.srt.bam $ILEXT | grep -v Mu | grep -v TE | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} 
( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${TE}.W_count.tsv
done
touch ${SAMP}_${TE}.W_count.tsv
done
done

# Check if W count results are intact
for TE in TE{1..32} ; do
for SAMP in $( cat ALL0510.txt ); do
a=$( wc -l ${SAMP}_${TE}_IL.txt | cut -d " " -f 1 )
b=$( wc -l ${SAMP}_${TE}.W_count.tsv | cut -d " " -f 1 )
c=$( wc -l ${SAMP}_${TE}.MASK_W_count.tsv | cut -d " " -f 1 )
if [[ ! a -eq b ]]; then
rm -f ${SAMP}_${TE}.W_count.tsv
for IL in $( cat ${SAMP}_${TE}_IL.txt ); do
ILEXT=$( echo $IL | awk '{match($0,/(chr[1-9AB]+):([0-9]+)-([0-9]+)$/,a); print a[1]":"a[2]-1001"-"a[3]+1000}' )
W=$( samtools view -q 1 -f 66 ${SAMP}.TE_05082022.srt.bam $ILEXT | grep -v Mu | grep -v TE | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} 
( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${TE}.W_count.tsv
done
fi
if [[ ! a -eq c ]]; then
cat ${SAMP}_${TE}_A_IL_BMASKED.bed | \
awk '$2>=1000 {print $1"\t"$2-1000"\t"$3+1000} $2<1000 {print $1"\t0\t"$3+1000}' | \
bedtools sort -i - > ${SAMP}_${TE}_A_IL_BMASKED.region
samtools view -h -@ 8 -L ${SAMP}_${TE}_A_IL_BMASKED.region ${SAMP}.TE_05082022.srt.bam | \
samtools fastq -1 ${SAMP}_${TE}.TMP_MASK_1.fq -2 ${SAMP}_${TE}.TMP_MASK_2.fq
grep "@" ${SAMP}_${TE}.TMP_MASK_1.fq | sort > ${SAMP}_${TE}.TMP_MASK_1.id
grep "@" ${SAMP}_${TE}.TMP_MASK_2.fq | sort > ${SAMP}_${TE}.TMP_MASK_2.id
comm -1 -2 ${SAMP}_${TE}.TMP_MASK_1.id ${SAMP}_${TE}.TMP_MASK_2.id > ${SAMP}_${TE}.TMP_MASK_12.id
sed -i 's/@//g' ${SAMP}_${TE}.TMP_MASK_12.id
seqtk subseq ${SAMP}_${TE}.TMP_MASK_1.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq
seqtk subseq ${SAMP}_${TE}.TMP_MASK_2.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq
bwa mem -t 7 DVSA_TE.05122022.fasta ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq \
${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq | samtools sort -o ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
samtools index ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
rm -f ${SAMP}_${TE}.MASK_W_count.tsv
for IL in $( cat ${SAMP}_${TE}_IL.txt ); do
if [[ ${IL} == *"A"* ]]; then
W=$( samtools view -q 1 -f 66 ${SAMP}_${TE}.PAIRED.TMP_MASK.bam $IL | \
grep -v TE | grep -v Mu | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); 
match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} ( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${TE}.MASK_W_count.tsv
fi
done
rm -f ${SAMP}_${TE}.TMP_MASK_1.fq 
rm -f ${SAMP}_${TE}.TMP_MASK_2.fq
rm -f ${SAMP}_${TE}.PAIRED.TMP_MASK.bam*
cat ${SAMP}_${TE}_B_IL_AMASKED.bed | \
awk '$2>=1000 {print $1"\t"$2-1000"\t"$3+1000} $2<1000 {print $1"\t0\t"$3+1000}' | \
bedtools sort -i - > ${SAMP}_${TE}_B_IL_AMASKED.region
samtools view -h -@ 8 -L ${SAMP}_${TE}_B_IL_AMASKED.region ${SAMP}.TE_05082022.srt.bam | \
samtools fastq -1 ${SAMP}_${TE}.TMP_MASK_1.fq -2 ${SAMP}_${TE}.TMP_MASK_2.fq
grep "@" ${SAMP}_${TE}.TMP_MASK_1.fq | sort > ${SAMP}_${TE}.TMP_MASK_1.id
grep "@" ${SAMP}_${TE}.TMP_MASK_2.fq | sort > ${SAMP}_${TE}.TMP_MASK_2.id
comm -1 -2 ${SAMP}_${TE}.TMP_MASK_1.id ${SAMP}_${TE}.TMP_MASK_2.id > ${SAMP}_${TE}.TMP_MASK_12.id
sed -i 's/@//g' ${SAMP}_${TE}.TMP_MASK_12.id
seqtk subseq ${SAMP}_${TE}.TMP_MASK_1.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq
seqtk subseq ${SAMP}_${TE}.TMP_MASK_2.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq
bwa mem -t 7 DVSB_TE.05122022.fasta ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq \
${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq | samtools sort -o ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
samtools index ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
for IL in $( cat ${SAMP}_${TE}_IL.txt ); do
if [[ ${IL} == *"B"* ]]; then
W=$( samtools view -q 1 -f 66 ${SAMP}_${TE}.PAIRED.TMP_MASK.bam $IL | \
grep -v TE | grep -v Mu | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); 
match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} ( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${TE}.MASK_W_count.tsv
fi
done
fi
done
done
```

 # Deal with reads mapped to partially masked regions

```bash
# Add masked intact and non-intact MULE regions in DVS_MU.bed
# Sort the DVS_TE.bed file
for N in {1..32} ; do
bedtools sort -i DVS_TE${N}.bed > DVS_TE${N}.srt.bed
done

# Remove DVS ILs from all the scanned ILs
for N in {1..32} ; do
bedtools intersect -v -wa -a ALL_TE${N}.2SPLIT_SUP.bed \
-b DVS_TE${N}.srt.bed > ALL_TE${N}.2SPLIT_SUP.nonDVS.bed
done

# Output surrounding 1000 bp sequences of filtered ILs and separate those on DVS_A and DVS_B
for N in {1..32} ; do
awk '$1~/A/ {print $1"\t"$2-500"\t"$3+500}' ALL_TE${N}.2SPLIT_SUP.nonDVS.bed \
> ALL_TE${N}.2SPLIT_SUP.nonDVS.chrA.SURR500.bed
awk '$1~/B/ {print $1"\t"$2-500"\t"$3+500}' ALL_TE${N}.2SPLIT_SUP.nonDVS.bed \
> ALL_TE${N}.2SPLIT_SUP.nonDVS.chrB.SURR500.bed
done

# Output surrounding 1000 bp sequences of filtered ILs and separate those on DVS_A and DVS_B
for N in {1..32} ; do
awk '{print $1":"$2"-"$3"\t"$1":"$2-500"-"$3+500}' ALL_TE${N}.2SPLIT_SUP.nonDVS.bed \
> TE${N}_IL_ILSURROUNDING.tsv
done

# Map the surrounding sequences of DVS_A ILs to unmasked DVS_B to check if their allelic regions are 
# masked in the modified DVS reference genome
for N in {1..32} ; do
bedtools getfasta -fi /zfs/socbd/bwu4/YU/annotation/CK2021.60.corrected.fasta \
-fo ALL_TE${N}.2SPLIT_SUP.nonDVS.chrA.SURR500.fasta \
-bed ALL_TE${N}.2SPLIT_SUP.nonDVS.chrA.SURR500.bed
bedtools getfasta -fi /zfs/socbd/bwu4/YU/annotation/CK2021.60.corrected.fasta \
-fo ALL_TE${N}.2SPLIT_SUP.nonDVS.chrB.SURR500.fasta \
-bed ALL_TE${N}.2SPLIT_SUP.nonDVS.chrB.SURR500.bed
minimap2 -x asm20 -t 8 CK2021.60.corrected.B.fasta \
ALL_TE${N}.2SPLIT_SUP.nonDVS.chrA.SURR500.fasta > ALL_TE${N}.chrA_IL_chrB_SURR500.paf
minimap2 -x asm20 -t 8 CK2021.60.corrected.A.fasta \
ALL_TE${N}.2SPLIT_SUP.nonDVS.chrB.SURR500.fasta > ALL_TE${N}.chrB_IL_chrA_SURR500.paf
cat ALL_TE${N}.chrA_IL_chrB_SURR500.paf | \
awk '$12>0 && $4-$3 > 200 {print $6"\t"$8"\t"$9"\t"$1}' > ALL_TE${N}.chrA_IL_SURR500_DVSB.bed
cat ALL_TE${N}.chrB_IL_chrA_SURR500.paf | \
awk '$12>0 && $4-$3 > 200 {print $6"\t"$8"\t"$9"\t"$1}' > ALL_TE${N}.chrB_IL_SURR500_DVSA.bed
bedtools intersect -wo -a ALL_TE${N}.chrA_IL_SURR500_DVSB.bed \
-b CK2021.0512_TE32_MU12_MASKED.bed > ALL_TE${N}.chrA_IL_SURR500_DVSB_MASKED.intersect.txt
bedtools intersect -wo -a ALL_TE${N}.chrB_IL_SURR500_DVSA.bed \
-b CK2021.0512_TE32_MU12_MASKED.bed > ALL_TE${N}.chrB_IL_SURR500_DVSA_MASKED.intersect.txt
done

for TE in TE{1..32}; do
cut -f 1,2,3,4 ALL_${TE}.chrA_IL_SURR500_DVSB_MASKED.intersect.txt | \
sort -u > ${TE}.A_IL_BMASKED.tsv
cut -f 1,2,3,4 ALL_${TE}.chrB_IL_SURR500_DVSA_MASKED.intersect.txt | \
sort -u > ${TE}.B_IL_AMASKED.tsv
done

# Prepare references including DVSA or DVSB and the modified transposon sequences
cdbfasta CK2021.05072022.MASKED.fasta
grep ">" CK2021.05072022.MASKED.fasta | grep -v B | sed "s/>//g" > DVSA_TE.05122022.gi
cat DVSA_TE.05122022.gi | cdbyank CK2021.05072022.MASKED.fasta.cidx - > DVSA_TE.05122022.fasta
bwa index DVSA_TE.05122022.fasta
grep ">" CK2021.05072022.MASKED.fasta | grep -v A | sed "s/>//g" > DVSB_TE.05122022.gi
cat DVSB_TE.05122022.gi | cdbyank CK2021.05072022.MASKED.fasta.cidx - > DVSB_TE.05122022.fasta
bwa index DVSB_TE.05122022.fasta

for TE in TE{1..32}; do
for SAMP in $( cat ALL0510.txt ); do
./FETCH_SAMP_ILS.py ${TE} $SAMP
cat ${SAMP}_${TE}_A_IL_BMASKED.bed | \
awk '$2>=1000 {print $1"\t"$2-1000"\t"$3+1000} $2<1000 {print $1"\t0\t"$3+1000}' | \
bedtools sort -i - > ${SAMP}_${TE}_A_IL_BMASKED.region
samtools view -h -@ 8 -L ${SAMP}_${TE}_A_IL_BMASKED.region ${SAMP}.TE_05082022.srt.bam | \
samtools fastq -1 ${SAMP}_${TE}.TMP_MASK_1.fq -2 ${SAMP}_${TE}.TMP_MASK_2.fq
grep "@" ${SAMP}_${TE}.TMP_MASK_1.fq | sort > ${SAMP}_${TE}.TMP_MASK_1.id
grep "@" ${SAMP}_${TE}.TMP_MASK_2.fq | sort > ${SAMP}_${TE}.TMP_MASK_2.id
comm -1 -2 ${SAMP}_${TE}.TMP_MASK_1.id ${SAMP}_${TE}.TMP_MASK_2.id > ${SAMP}_${TE}.TMP_MASK_12.id
sed -i 's/@//g' ${SAMP}_${TE}.TMP_MASK_12.id
seqtk subseq ${SAMP}_${TE}.TMP_MASK_1.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq
seqtk subseq ${SAMP}_${TE}.TMP_MASK_2.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq
bwa mem -t 7 DVSA_TE.05122022.fasta ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq \
${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq | samtools sort -o ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
samtools index ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
rm -f ${SAMP}_${TE}.MASK_W_count.tsv
for IL in $( cat ${SAMP}_${TE}_IL.txt ); do
if [[ ${IL} == *"A"* ]]; then
W=$( samtools view -q 1 -f 66 ${SAMP}_${TE}.PAIRED.TMP_MASK.bam $IL | \
grep -v TE | grep -v Mu | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); 
match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} ( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${TE}.MASK_W_count.tsv
fi
done
rm -f ${SAMP}_${TE}.TMP_MASK_1.fq 
rm -f ${SAMP}_${TE}.TMP_MASK_2.fq
rm -f ${SAMP}_${TE}.PAIRED.TMP_MASK.bam*
cat ${SAMP}_${TE}_B_IL_AMASKED.bed | \
awk '$2>=1000 {print $1"\t"$2-1000"\t"$3+1000} $2<1000 {print $1"\t0\t"$3+1000}' | \
bedtools sort -i - > ${SAMP}_${TE}_B_IL_AMASKED.region
samtools view -h -@ 8 -L ${SAMP}_${TE}_B_IL_AMASKED.region ${SAMP}.TE_05082022.srt.bam | \
samtools fastq -1 ${SAMP}_${TE}.TMP_MASK_1.fq -2 ${SAMP}_${TE}.TMP_MASK_2.fq
grep "@" ${SAMP}_${TE}.TMP_MASK_1.fq | sort > ${SAMP}_${TE}.TMP_MASK_1.id
grep "@" ${SAMP}_${TE}.TMP_MASK_2.fq | sort > ${SAMP}_${TE}.TMP_MASK_2.id
comm -1 -2 ${SAMP}_${TE}.TMP_MASK_1.id ${SAMP}_${TE}.TMP_MASK_2.id > ${SAMP}_${TE}.TMP_MASK_12.id
sed -i 's/@//g' ${SAMP}_${TE}.TMP_MASK_12.id
seqtk subseq ${SAMP}_${TE}.TMP_MASK_1.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq
seqtk subseq ${SAMP}_${TE}.TMP_MASK_2.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq
bwa mem -t 7 DVSB_TE.05122022.fasta ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq \
${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq | samtools sort -o ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
samtools index ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
for IL in $( cat ${SAMP}_${TE}_IL.txt ); do
if [[ ${IL} == *"B"* ]]; then
W=$( samtools view -q 1 -f 66 ${SAMP}_${TE}.PAIRED.TMP_MASK.bam $IL | \
grep -v TE | grep -v Mu | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); 
match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} ( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${TE}.MASK_W_count.tsv
fi
done
done
done
```

### # FETCH_SAMP_ILS.py content

```python
#!/usr/bin/env python3

import pandas as pd
import sys

TE = sys.argv[1] # TE name
SAMP = sys.argv[2] # Sample name in the genotype table
GTDF = pd.read_table(f'{TE}_GENOTYPE.tsv', header=0)
GTDF.loc[GTDF[SAMP]==1,'name'].to_csv(f'{SAMP}_{TE}_IL.txt',header=False,index=False)
ID_CONVER = pd.read_table(f'{TE}_IL_ILSURROUNDING.tsv',names=['name','surr'])

A_IL = pd.read_table(f'{TE}.A_IL_BMASKED.tsv',names=['chr','Left','Right','surr'])
A_IL = A_IL.merge(ID_CONVER,on='surr')
A_IL.loc[A_IL['name'].isin(GTDF.loc[GTDF[SAMP]==1,'name']),['chr','Left','Right']].to_csv(f'{SAMP}_{TE}_A_IL_BMASKED.bed',sep="\t",header=False,index=False)
A_IL.loc[A_IL['name'].isin(GTDF.loc[GTDF[SAMP]==1,'name']),['surr']].to_csv(f'{SAMP}_{TE}_A_IL_SURR.bed',header=False,index=False)

B_IL = pd.read_table(f'{TE}.B_IL_AMASKED.tsv',names=['chr','Left','Right','surr'])
ID_CONVER = pd.read_table(f'{TE}_IL_ILSURROUNDING.tsv',names=['name','surr'])
B_IL = B_IL.merge(ID_CONVER,on='surr')
B_IL.loc[B_IL['name'].isin(GTDF.loc[GTDF[SAMP]==1,'name']),['chr','Left','Right']].to_csv(f'{SAMP}_{TE}_B_IL_AMASKED.bed',sep="\t",header=False,index=False)
B_IL.loc[B_IL['name'].isin(GTDF.loc[GTDF[SAMP]==1,'name']),['surr']].to_csv(f'{SAMP}_{TE}_B_IL_SURR.bed',header=False,index=False)
```

 ### Count wild-type reads mapped to the partially masked allelic regions  for PRECISE_ILs not detected by the old method

```bash
# Screen for the Precise ILs without wild-type read count statistics in previsous analysis
for TE in TE{1..32}; do
cut -f 1,2 GENOTYPE/${TE}_0529_GENOTYPE.tsv | awk 'NR>1 && $2~/DVS/' | \
sed "s/:/\t/g;s/-/\t/g" | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' \
> GENOTYPE/${TE}_0529_GENOTYPE_NONDVS.bed

awk 'NR>1 {print $230}' GENOTYPE/${TE}_GENOTYPE.tsv | \
sed "s/:/\t/g;s/-/\t/g" | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' \
> GENOTYPE/${TE}_GENOTYPE.bed

bedtools intersect -v -a GENOTYPE/${TE}_0529_GENOTYPE_NONDVS.bed \
-b GENOTYPE/${TE}_GENOTYPE.bed > TEMP/${TE}_0529_GENOTYPE_NONDVS_LEFT_TMP.bed 
done

# Output surrounding 1000 bp sequences of rest non-DVS ILs and separate those on DVS_A and DVS_B
for N in {1..32} ; do
awk '$1~/A/ {print $1"\t"$2-500"\t"$3+500}' TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.bed \
> TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.chrA.SURR500.bed
awk '$1~/B/ {print $1"\t"$2-500"\t"$3+500}' TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.bed \
> TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.chrB.SURR500.bed
done

# Output surrounding 1000 bp sequences of filtered ILs and separate those on DVS_A and DVS_B
for N in {1..32} ; do
awk '{print $1":"$2"-"$3"\t"$1":"$2-500"-"$3+500}' TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.bed \
> TEMP/TE${N}_IL_ILSURROUNDING_TMP.tsv
done

# Map the surrounding sequences of DVS_A ILs to unmasked DVS_B to check if their allelic regions are 
# masked in the modified DVS reference genome
for N in {1..32} ; do
bedtools getfasta -fi /zfs/socbd/bwu4/YU/annotation/CK2021.60.corrected.fasta \
-fo TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.chrA.SURR500.fasta \
-bed TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.chrA.SURR500.bed
bedtools getfasta -fi /zfs/socbd/bwu4/YU/annotation/CK2021.60.corrected.fasta \
-fo TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.chrB.SURR500.fasta \
-bed TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.chrB.SURR500.bed
minimap2 -x asm20 -t 8 CK2021.60.corrected.B.fasta \
TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.chrA.SURR500.fasta > TEMP/TE${N}_REST.chrA_IL_chrB_SURR500_TMP.paf
minimap2 -x asm20 -t 8 CK2021.60.corrected.A.fasta \
TEMP/TE${N}_0529_GENOTYPE_NONDVS_LEFT_TMP.chrB.SURR500.fasta > TEMP/TE${N}_REST.chrB_IL_chrA_SURR500_TMP.paf
cat TEMP/TE${N}_REST.chrA_IL_chrB_SURR500_TMP.paf | \
awk '$12>0 && $4-$3 > 200 {print $6"\t"$8"\t"$9"\t"$1}' > TEMP/TE${N}_REST.chrA_IL_SURR500_DVSB_TMP.bed
cat TEMP/TE${N}_REST.chrB_IL_chrA_SURR500_TMP.paf | \
awk '$12>0 && $4-$3 > 200 {print $6"\t"$8"\t"$9"\t"$1}' > TEMP/TE${N}_REST.chrB_IL_SURR500_DVSA_TMP.bed
bedtools intersect -wo -a TEMP/TE${N}_REST.chrA_IL_SURR500_DVSB_TMP.bed \
-b CK2021.0512_TE32_MU12_MASKED.bed > TEMP/TE${N}_REST.chrA_IL_SURR500_DVSB_MASKED.intersect.txt
bedtools intersect -wo -a TEMP/TE${N}_REST.chrB_IL_SURR500_DVSA_TMP.bed \
-b CK2021.0512_TE32_MU12_MASKED.bed > TEMP/TE${N}_REST.chrB_IL_SURR500_DVSA_MASKED.intersect.txt
done

for N in {1..32}; do
cut -f 1,2,3,4 TEMP/TE${N}_REST.chrA_IL_SURR500_DVSB_MASKED.intersect.txt | \
sort -u > TE${N}_REST.A_IL_BMASKED.tsv
cut -f 1,2,3,4 TEMP/TE${N}_REST.chrB_IL_SURR500_DVSA_MASKED.intersect.txt | \
sort -u > TE${N}_REST.B_IL_AMASKED.tsv
done

for SAMP in $( cat ALL0510.txt ); do
for TE in TE{1..32} ;do
./NEW_FETCH_SAMP_ILS.py $TE $SAMP
done
deon

for SAMP in $( cat ALL0510.txt ); do
sed "s/SAMP/${SAMP}/g" TEMPLATE_MASK_W_COUNT.pbs > ${SAMP}_MASK_W_COUNT.pbs
qsub ${SAMP}_MASK_W_COUNT.pbs
done

```

### # Content of MASK_W_COUNT.sh

```bash
#!/bin/bash

SAMP=$1
THREADS=7
for TE in TE{1..32}; do
cat TEMP/A_IL_BMASKED/${SAMP}_${TE}_REST_A_IL_BMASKED.bed | \
awk '$2>=500 {print $1"\t"$2-500"\t"$3+500} $2<500 {print $1"\t0\t"$3+500}' | \
bedtools sort -i - > ${SAMP}_${TE}.REST_A_IL_BMASKED_TMP.region
samtools view -h -q 1 -@ $THREADS -L ${SAMP}_${TE}.REST_A_IL_BMASKED_TMP.region BAM/${SAMP}.TE_05082022.srt.bam | \
samtools fastq -1 ${SAMP}_${TE}.TMP_MASK_1.fq -2 ${SAMP}_${TE}.TMP_MASK_2.fq
grep "@" ${SAMP}_${TE}.TMP_MASK_1.fq | sort > ${SAMP}_${TE}.TMP_MASK_1.id
grep "@" ${SAMP}_${TE}.TMP_MASK_2.fq | sort > ${SAMP}_${TE}.TMP_MASK_2.id
comm -1 -2 ${SAMP}_${TE}.TMP_MASK_1.id ${SAMP}_${TE}.TMP_MASK_2.id > ${SAMP}_${TE}.TMP_MASK_12.id
sed -i 's/@//g' ${SAMP}_${TE}.TMP_MASK_12.id
seqtk subseq ${SAMP}_${TE}.TMP_MASK_1.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq
seqtk subseq ${SAMP}_${TE}.TMP_MASK_2.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq
bwa mem -t ${THREADS} DVSA_TE.05122022.fasta ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq \
${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq | samtools sort -o ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
samtools index -@ $THREADS ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
rm -f ${SAMP}_${TE}_REST.MASK_W_count.tsv
touch ${SAMP}_${TE}_REST.MASK_W_count.tsv
for IL in $( cat TEMP/${SAMP}_${TE}_REST_NONDVS_IL.txt ); do
if [[ ${IL} == *"A"* ]]; then
W=$( samtools view -q 1 -f 66 ${SAMP}_${TE}.PAIRED.TMP_MASK.bam $IL | \
grep -v TE | grep -v Mu | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); 
match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} ( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${TE}_REST.MASK_W_count.tsv
fi
done
rm -f ${SAMP}_${TE}.*TMP* 
cat TEMP/B_IL_AMASKED/${SAMP}_${TE}_REST_B_IL_AMASKED.bed | \
awk '$2>=500 {print $1"\t"$2-500"\t"$3+1000} $2<500 {print $1"\t0\t"$3+500}' | \
bedtools sort -i - > ${SAMP}_${TE}.B_IL_AMASKED_TMP.region
samtools view -h -q 1 -@ ${THREADS} -L ${SAMP}_${TE}.B_IL_AMASKED_TMP.region BAM/${SAMP}.TE_05082022.srt.bam | \
samtools fastq -1 ${SAMP}_${TE}.TMP_MASK_1.fq -2 ${SAMP}_${TE}.TMP_MASK_2.fq
grep "@" ${SAMP}_${TE}.TMP_MASK_1.fq | sort > ${SAMP}_${TE}.TMP_MASK_1.id
grep "@" ${SAMP}_${TE}.TMP_MASK_2.fq | sort > ${SAMP}_${TE}.TMP_MASK_2.id
comm -1 -2 ${SAMP}_${TE}.TMP_MASK_1.id ${SAMP}_${TE}.TMP_MASK_2.id > ${SAMP}_${TE}.TMP_MASK_12.id
sed -i 's/@//g' ${SAMP}_${TE}.TMP_MASK_12.id
seqtk subseq ${SAMP}_${TE}.TMP_MASK_1.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq
seqtk subseq ${SAMP}_${TE}.TMP_MASK_2.fq ${SAMP}_${TE}.TMP_MASK_12.id | paste - - - - | \
sort -k1,1 -S 3G | tr '\t' '\n' > ${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq
bwa mem -t ${THREADS} DVSB_TE.05122022.fasta ${SAMP}_${TE}.PAIRED.TMP_MASK_1.fq \
${SAMP}_${TE}.PAIRED.TMP_MASK_2.fq | samtools sort -o ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
samtools index -@ ${THREADS} ${SAMP}_${TE}.PAIRED.TMP_MASK.bam
for IL in $( cat TEMP/${SAMP}_${TE}_REST_NONDVS_IL.txt ); do
if [[ ${IL} == *"B"* ]]; then
W=$( samtools view -q 1 -f 66 ${SAMP}_${TE}.PAIRED.TMP_MASK.bam $IL | \
grep -v TE | grep -v Mu | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); 
match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} ( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${TE}_REST.MASK_W_count.tsv
fi
done
rm -f ${SAMP}_${TE}.*TMP* 
done

```

### # Content of W_COUNT.sh

```bash
#!/bin/bash
SAMP=$1
THREADS=$2
for TE in TE{1..32} ; do
rm -f ${SAMP}_${TE}_REST.W_count.tsv
touch ${SAMP}_${TE}_REST.W_count.tsv
for IL in $( cat TEMP/${SAMP}_${TE}_REST_NONDVS_IL.txt ); do
ILEXT=$( echo $IL | awk '{match($0,/(chr[1-9AB]+):([0-9]+)-([0-9]+)$/,a); print a[1]":"a[2]-1001"-"a[3]+1000}' )
W=$( samtools view -@ ${THREADS} -q 1 -f 66 BAM/${SAMP}.TE_05082022.srt.bam $ILEXT | grep -v Mu | grep -v TE | \
awk -v x=$IL '{match(x,/:([0-9]+)-([0-9]+)/,a); match($6,/([0-9]+)M/,b); match($0,/MC:Z:([0-9]+[SH])?([0-9]+)M/,c)} 
( $4<=a[1]-40 && $8+c[2]>=a[2]+40) || ($8<=a[1]-40 && $4+b[1]>=a[2]+40) {print $0}' | wc -l )
echo -e ${IL}"\t"${W} >> ${SAMP}_${TE}_REST.W_count.tsv
done
done

```

### # Content of TEMPLATE_MASK_W_COUNT.pbs

```bash
#PBS -N SAMP
#PBS -l select=1:ncpus=8:mem=30gb:interconnect=1g,walltime=168:00:00
#PBS -j oe

source ~/.bashrc
cd /scratch1/bwu4/TE_SCAN
./MASK_W_COUNT.sh SAMP
```

### # NEW_FETCH_SAMP_ILS.py content

```python
#!/usr/bin/env python3

import pandas as pd
import sys

TE = sys.argv[1] # TE name
SAMP = sys.argv[2] # Sample name in the genotype table
GTDF = pd.read_table(f'GENOTYPE/{TE}_0529_GENOTYPE.tsv', header=0)
REST = pd.read_table(f'TEMP/{TE}_0529_GENOTYPE_NONDVS_LEFT_TMP.bed',names=['CHR','START','END','name'], header=None)
REST = REST.loc[REST['name'].str.contains('chr')].copy()
GTDF = GTDF.loc[(GTDF[SAMP]==1) & (GTDF['name'].isin(REST['name']))]
GTDF['name'].to_csv(f'TEMP/{SAMP}_{TE}_REST_NONDVS_IL.txt',header=False,index=False)
ID_CONVER = pd.read_table(f'TEMP/{TE}_IL_ILSURROUNDING_TMP.tsv',names=['name','surr'])

A_IL = pd.read_table(f'{TE}_REST.A_IL_BMASKED.tsv',names=['chr','Left','Right','surr'])
A_IL = A_IL.merge(ID_CONVER,on='surr')
A_IL.loc[A_IL['name'].isin(GTDF['name']),['chr','Left','Right']].to_csv(f'TEMP/{SAMP}_{TE}_REST_A_IL_BMASKED.bed',sep="\t",header=False,index=False)
A_IL.loc[A_IL['name'].isin(GTDF['name']),['surr']].to_csv(f'TEMP/{SAMP}_{TE}_REST_A_IL_SURR.bed',header=False,index=False)

B_IL = pd.read_table(f'{TE}_REST.B_IL_AMASKED.tsv',names=['chr','Left','Right','surr'])
B_IL = B_IL.merge(ID_CONVER,on='surr')
B_IL.loc[B_IL['name'].isin(GTDF['name']),['chr','Left','Right']].to_csv(f'TEMP/{SAMP}_{TE}_REST_B_IL_AMASKED.bed',sep="\t",header=False,index=False)
B_IL.loc[B_IL['name'].isin(GTDF['name']),['surr']].to_csv(f'TEMP/{SAMP}_{TE}_REST_B_IL_SURR.bed',header=False,index=False)
```

## # Statistics on IL counts on TE1-32 for SWO, PUM, and MAN

```python
#!/usr/bin/env python3

import pandas as pd

LIST=[]
for i in range(1,33):
    LIST.append("TE"+str(i))

ID = pd.read_table('SRR_SPECIES.tsv', header=0)
SWO = ID.loc[ID['Species']=="SWO",'ID'].to_list()
MAN = ID.loc[ID['Species']=="MAN",'ID'].to_list()
PUM = ID.loc[ID['Species']=="PUM",'ID'].to_list()

df = pd.DataFrame()
df['TE'] = LIST

for TE in LIST:
    GENOTYPE = pd.read_table(TE+'_GENOTYPE.tsv', header=0)
    for SPE in ['SWO', 'MAN', 'PUM']:
        df.loc[df['TE']==TE, SPE+'_mean'] = GENOTYPE[globals()[SPE]].sum(axis=0).mean()
        df.loc[df['TE']==TE, SPE+'_min'] = GENOTYPE[globals()[SPE]].sum(axis=0).min()
        df.loc[df['TE']==TE, SPE+'_max'] = GENOTYPE[globals()[SPE]].sum(axis=0).max()
        df.loc[df['TE']==TE, SPE+'_std'] = GENOTYPE[globals()[SPE]].sum(axis=0).std()
    GENOTYPE_T = GENOTYPE.set_index('name').T
    GENOTYPE_T.reset_index(level=0,inplace=True)
    GENOTYPE_T = GENOTYPE_T.rename(columns={"index": "ID"})
    GENOTYPE_T = GENOTYPE_T.merge(ID,on='ID')
    SUM = GENOTYPE_T.set_index('ID').groupby('Species').sum()
    SUMT = SUM.T
    df.loc[df['TE']==TE, 'MAN_UNIQ_IL'] = len(SUMT.loc[(SUMT['MAN']>0) & (SUMT['PUM']==0) & (SUMT['SWO']==0)].index)
    df.loc[df['TE']==TE, 'PUM_UNIQ_IL'] = len(SUMT.loc[(SUMT['PUM']>0) & (SUMT['MAN']==0) & (SUMT['SWO']==0)].index)
    df.loc[df['TE']==TE, 'SWO_UNIQ_IL'] = len(SUMT.loc[(SUMT['SWO']>0) & (SUMT['MAN']==0) & (SUMT['PUM']==0)].index)

df.to_csv('TE1-32_Species_IL_Statistics.tsv', sep="\t", header=True, index=False)
```

### TE classification

```bash
qsub -I -l select=1:ncpus=24:mem=123gb:ngpus=1:gpu_model=k40:interconnect=fdr,walltime=72:00:00
conda activate py36
cd /zfs/socbd/bwu4/YU/annotation/FOUR_TE/TE_IN_CITRUS_GENOMES/SWO
DeepTE/DeepTE.py -d working_dir -o DeepTE_output -i ALLTE_REP.fasta -sp P -m P
```

### # Correct DVS ILs located on the allelic chromosome

```bash
for N in {1..32}; do
# Output upstream and downstream 100 bp regions of the scanned ILs
awk '$1~/:/ {match($1, /(chr[1-9AB]*):([0-9]*)-([0-9]*)/, a); print a[1]"\t"a[3]-101"\t"a[3]"\t"$1}' \
/scratch1/bwu4/TE_SCAN/TE${N}_NEW_GENOTYPE.tsv > TE${N}_IL_LEFT100.bed
bedtools getfasta -nameOnly -fi CK2021.60.corrected.fasta -fo TE${N}_IL_LEFT100.fasta -bed TE${N}_IL_LEFT100.bed
awk '$1~/:/ {match($1, /(chr[1-9AB]*):([0-9]*)-([0-9]*)/, a); print a[1]"\t"a[2]-1"\t"a[2]+100"\t"$1}' \
/scratch1/bwu4/TE_SCAN/TE${N}_NEW_GENOTYPE.tsv > TE${N}_IL_RIGHT100.bed
bedtools getfasta -nameOnly -fi CK2021.60.corrected.fasta -fo TE${N}_IL_RIGHT100.fasta -bed TE${N}_IL_RIGHT100.bed
makeblastdb -in TE${N}_IL_LEFT100.fasta -parse_seqids -out TE${N}_IL_LEFT100 -dbtype nucl
makeblastdb -in TE${N}_IL_RIGHT100.fasta -parse_seqids -out TE${N}_IL_RIGHT100 -dbtype nucl

# Output upstream and downstream 100 bp regions of the DVS intact and unintact members
awk '{print $1"\t"$2-101"\t"$3+100}' NEW${N}_UNINTACT_COPIES.bed > NEW${N}_UNINTACT_EXT100.bed
bedtools getfasta -fi CK2021.60.corrected.fasta -fo NEW${N}_UNINTACT_EXT100.fasta \
-bed NEW${N}_UNINTACT_EXT100.bed
awk '{print $1"\t"$2-101"\t"$2}' NEW${N}_INTACT_COPIES.bed > NEW${N}_INTACT_LEFT100.bed
awk '{print $1"\t"$3-1"\t"$3+100}' NEW${N}_INTACT_COPIES.bed > NEW${N}_INTACT_RIGHT100.bed
bedtools getfasta -fi CK2021.60.corrected.fasta -fo NEW${N}_INTACT_L100.fasta \
-bed NEW${N}_INTACT_LEFT100.bed
bedtools getfasta -fi CK2021.60.corrected.fasta -bed NEW${N}_INTACT_RIGHT100.bed > NEW${N}_INTACT_R100.fasta

blastn -query NEW${N}_INTACT_L100.fasta -db TE${N}_IL_RIGHT100 -outfmt 6 | awk '($1~/A/ && $2~/B/) || ($1~/B/ && $2~/A/)' | awk '$8>85 && $3>95 && $10<15 {print $1"\t"$2}' > TE${N}_AB_FILT_IL.tsv
blastn -query NEW${N}_INTACT_R100.fasta -db TE${N}_IL_RIGHT100 -outfmt 6 | awk '($1~/A/ && $2~/B/) || ($1~/B/ && $2~/A/)' | awk '$7<15 && $3>95 && $9<15  {print $1"\t"$2}' >> TE${N}_AB_FILT_IL.tsv
blastn -query NEW${N}_INTACT_L100.fasta -db TE${N}_IL_LEFT100 -outfmt 6 | awk '($1~/A/ && $2~/B/) || ($1~/B/ && $2~/A/)' | awk '$8>85 && $3>95 && $10>85  {print $1"\t"$2}' >> TE${N}_AB_FILT_IL.tsv
blastn -query NEW${N}_INTACT_R100.fasta -db TE${N}_IL_LEFT100 -outfmt 6 | awk '($1~/A/ && $2~/B/) || ($1~/B/ && $2~/A/)' | awk '$7<15 && $3>95 && $9>85 {print $1"\t"$2}' >> TE${N}_AB_FILT_IL.tsv
done

for N in {1..32}; do
sort -u TE${N}_AB_FILT_IL.tsv > TE${N}_AB_FILT_IL.tmp
mv TE${N}_AB_FILT_IL.tmp TE${N}_AB_FILT_IL.tsv
done 

for TE in TE{1..32}; do ./AB_FILT_IL.py $TE ; done
```

### ###### Content of AB_FILT_IL.py

```python
#!/usr/bin/env python3

import pandas as pd
import sys

TE=sys.argv[1]
DVS_IL = pd.read_table(f'{TE}_AB_FILT_IL.tsv', names=['DVS','name'])
GENOTYPE = pd.read_table(f'{TE}_NEW_GENOTYPE.tsv', header=0)
DVS_IL = DVS_IL.merge(GENOTYPE[['name', 'DVS']],on='name')

GENOTYPE.loc[~GENOTYPE['name'].isin(DVS_IL.loc[DVS_IL['DVS_y'].str.contains('NOT'),'name'])].to_csv(f'{TE}_GENOTYPE_FILT.tsv',header=True,sep="\t")
```

### ###### Calculate Mutant cell ratio

```bash
for TE in TE{1..32}; do
./Calculate_MRATIO.py $TE
done
```

### ###### Content of Calculate_MRATIO.py

```python
#!/usr/bin/env python3

import pandas as pd
import sys

TE=sys.argv[1]
W_COUNT = pd.read_table(f'{TE}_GENOTYPE_FILT.tsv', header=0)
MASKW_COUNT = pd.read_table(f'{TE}_GENOTYPE_FILT.tsv', header=0)
SAMP = pd.read_table('ALL0510.txt', header=None)

for SRR in SAMP[0]:
    W_COUNT[SRR] = 0
    W = pd.read_table(f'{SRR}_{TE}.W_count.tsv', names=['name','W_COUNT'])
    for IL in W['name']:
        if IL in W_COUNT['name'].values:
            W_COUNT.loc[W_COUNT['name']==IL,SRR] = W.loc[W['name']==IL,'W_COUNT'].to_list()[0]

for SRR in SAMP[0]:
    MASKW_COUNT[SRR] = 0
    MASKW = pd.read_table(f'{SRR}_{TE}.MASK_W_count.tsv', names=['name','W_COUNT'])
    for IL in MASKW['name']:
        if IL in MASKW_COUNT['name'].values:
            MASKW_COUNT.loc[MASKW_COUNT['name']==IL,SRR] = MASKW.loc[MASKW['name']==IL,'W_COUNT'].to_list()[0]

MASKA = pd.read_table(f'{TE}.A_IL_BMASKED.tsv',header=None)[3].drop_duplicates()
MASKB = pd.read_table(f'{TE}.B_IL_AMASKED.tsv',header=None)[3].drop_duplicates()
MASK = pd.concat([MASKA,MASKB], ignore_index=True)

W_COUNT = W_COUNT.set_index(['name','DVS'])
MASKW_COUNT = MASKW_COUNT.set_index(['name','DVS'])
WMERGE_COUNT = W_COUNT + MASKW_COUNT
WMERGE_COUNT = WMERGE_COUNT.reset_index()
WMERGE_COUNT.to_csv(f'{TE}_MERGED_W_READ_COUNT.tsv', sep="\t", header=True, index=False)
WMERGE_COUNT_MASK = WMERGE_COUNT.loc[WMERGE_COUNT['name'].isin(MASK)].copy()
WMERGE_COUNT_NOTMASK = WMERGE_COUNT.loc[~WMERGE_COUNT['name'].isin(MASK)].copy()
WMERGE_COUNT_MASK.set_index(['name', 'DVS'], inplace=True)
WMERGE_COUNT_NOTMASK.set_index(['name', 'DVS'], inplace=True)

M_COUNT = pd.read_table(f'ALL_SAMP.{TE}.SPLIT_DISC_DVSMERGED.read_count.tsv', header=0)
M_COUNT = M_COUNT.loc[M_COUNT['name'].isin(WMERGE_COUNT['name'])]
M_COUNT_MASK = M_COUNT.loc[M_COUNT['name'].isin(MASK)].copy()
M_COUNT_NOTMASK = M_COUNT.loc[~M_COUNT['name'].isin(MASK)].copy()
M_COUNT_MASK.set_index(['name', 'DVS'], inplace=True)
M_COUNT_NOTMASK.set_index(['name', 'DVS'], inplace=True)

NOTMASK_MRATIO = M_COUNT_NOTMASK/2/(M_COUNT_NOTMASK/2+WMERGE_COUNT_NOTMASK+0.001)
NOTMASK_MRATIO.where(NOTMASK_MRATIO<1,1,inplace=True)
NOTMASK_MRATIO.reset_index(inplace=True)
MASK_MRATIO = M_COUNT_MASK/(M_COUNT_MASK/2+WMERGE_COUNT_MASK+0.001)
MASK_MRATIO.where(MASK_MRATIO<1,1,inplace=True)
MASK_MRATIO.reset_index(inplace=True)

pd.concat([NOTMASK_MRATIO, MASK_MRATIO], ignore_index=True).to_csv(f'{TE}_MRATIO.tsv', sep="\t", header=True, index=False)
```

### ###### Merge putative false-positive allelic ILs

```bash
cat TE1_GENOTYPE_FILT.tsv | awk \
'$232~/NOT_DVS/{match($2,/(chr[1-9AB]*):([0-9]+)-([0-9]+)/, a); print a[1]"\t"a[2]-50"\t"a[3]+50"\t"$2}' \
> TE1_NOTDVS_50SUR.bed

bedtools getfasta -nameOnly -fi CK2021.60.corrected.fasta -bed TE1_NOTDVS_50SUR.bed -fo TE1_NOTDVS_50SUR.fasta

makeblastdb -in TE1_NOTDVS_50SUR.fasta -dbtype nucl -parse_seqids -out TE1_NOTDVS_50SUR
blastn -query TE1_NOTDVS_50SUR.fasta -db TE1_NOTDVS_50SUR -outfmt 6 -num_threads 24 | awk '$1!=$2'
```

### # Scanning precise ILs

```bash
bwa index ALLTE_REP.fasta
bwa mem -t 20 ALLTE_REP.fasta CK_1.clean.fastq.gz  CK_2.clean.fastq.gz | samtools view -h -q1 -F4 > CK_ALLREP.sam
samtools sort -@ 20 -o CK_ALLREP.srt.bam CK_ALLREP.sam
samtools view CK_ALLREP.srt.bam | awk '$3~/TE1_/' | awk ' 
{if ($5>0 && $4==1 && $6~/S/) {match($6, /([0-9]+)S([0-9]+)M/, a); 
if ( a[1]>15 && a[2]>=30) {print ">"$1"\t"substr($10,1,a[1])}}}' | \
sed "s/\t/\n/g" > CK_TE1.split_reads_L.fasta

for N in {1..32} ; do
infoseq -only -name -length NEW${N}_REP.fasta | \
awk '$1~/TE/ {print $1":"$2"-"$2}' > TE${N}_R_Terminal.txt 
done

for TERMINAL in $( cat TE1_R_Terminal.txt ); do
samtools view CK_TE1.srt.bam $TERMINAL | awk '
{if ($5>0 && $6~/S/) {match($6, /([0-9]+)M([0-9]+)S/, a);
if ( a[1]>=30 && a[2]>15 ) {print ">"$1"\t"substr($10,a[1]+1)}}
}' | sed "s/\t/\n/g" >> CK_TE1.split_reads_R.fasta
done

bwa mem CK2021.05072022.MASKED.fasta CK_TE1.split_reads_L.fasta | \
samtools view -q 1 -F 16 | \
awk '{match($6, /([0-9]+)M/, a); print $3"\t"$4+a[1]-1"\t"$4+a[1]}' > CK_TE1.PRECISE.SPLIT_IL.bed

bwa mem CK2021.05072022.MASKED.fasta CK_TE1.split_reads_L.fasta | \
samtools view -q 1 -f 16 | \
awk '{match($6, /([0-9]+)M/, a); print $3"\t"$4-1"\t"$4}' >> CK_TE1.PRECISE.SPLIT_IL.bed

bwa mem CK2021.05072022.MASKED.fasta CK_TE1.split_reads_R.fasta | \
samtools view -q 1 -f 16 | \
awk '{match($6, /([0-9]+)M/, a); print $3"\t"$4+a[1]-2"\t"$4+a[1]-1}' >> CK_TE1.PRECISE.SPLIT_IL.bed

bwa mem CK2021.05072022.MASKED.fasta CK_TE1.split_reads_R.fasta | \
samtools view -q 1 -F 16 | \
awk '{match($6, /([0-9]+)M/, a); print $3"\t"$4-1"\t"$4}' >> CK_TE1.PRECISE.SPLIT_IL.bed

sort -u CK_TE1.PRECISE.SPLIT_IL.bed | bedtools sort -i - | \
bedtools merge -d 20 -i - > CK_TE1.PRECISE.SPLIT_IL_MERGED.bed

bedtools intersect -c -a CK_TE1.PRECISE.SPLIT_IL_MERGED.bed -b CK_TE1.PRECISE.SPLIT_IL.bed
bedtools intersect -c -a CK_TE1.PRECISE.SPLIT_IL_MERGED.bed -b TE1/CK.TE1.DISC_SUP.bed
bedtools intersect -c -a CK_TE1.PRECISE.SPLIT_IL_MERGED.bed -b TE1/CK.TE1.SPLIT_SUP.bed

```

### # Content of MAP_TO_ALLTEREP.sh

```bash
#!/bin/bash
SAMP=$1
TE=$2
bwa mem -t 20 ALLTE_REP.fasta ${SAMP}_1.clean.fastq  ${SAMP}_2.clean.fastq | \
samtools view -h -q1 -F4 > ${SAMP}_ALLREP.sam
samtools sort -@ 20 -o ${SAMP}_ALLREP.srt.bam ${SAMP}_ALLREP.sam

```

```bash
TE=TE30
for SAMP in $( grep RR ALL0510.txt ); do
pat=${TE}"_"
#samtools index ${SAMP}_ALLREP.srt.bam
samtools view ${SAMP}_ALLREP.srt.bam | awk -v x="$pat" '$3~x' | awk ' 
{if ($5>0 && $4==1 && $6~/S/) {match($6, /([0-9]+)S([0-9]+)M/, a); 
if ( a[1]>15 && a[2]>=30) {print ">"$1"\t"substr($10,1,a[1])}}}' | \
sed "s/\t/\n/g" > ${SAMP}_${TE}.split_reads_L.fasta
rm -f ${SAMP}_${TE}.split_reads_R.fasta
for TERMINAL in $( cat ${TE}_R_Terminal.txt ); do
samtools view ${SAMP}_ALLREP.srt.bam $TERMINAL | awk '
{if ($5>0 && $6~/S/) {match($6, /([0-9]+)M([0-9]+)S/, a);
if ( a[1]>=30 && a[2]>15 ) {print ">"$1"\t"substr($10,a[1]+1)}}
}' | sed "s/\t/\n/g" >> ${SAMP}_${TE}.split_reads_R.fasta
done
done

TE=TE32
for TE in TE{1..30}; do
for SAMP in $( cat ALL0510.txt ); do
seqkit rename ${SAMP}_${TE}.split_reads_L.fasta > ${SAMP}_${TE}.split_reads_L.fasta1
mv ${SAMP}_${TE}.split_reads_L.fasta1 ${SAMP}_${TE}.split_reads_L.fas
seqkit rename ${SAMP}_${TE}.split_reads_R.fasta > ${SAMP}_${TE}.split_reads_R.fasta1
mv ${SAMP}_${TE}.split_reads_R.fasta1 ${SAMP}_${TE}.split_reads_R.fas
done
done
```

### # Content of OUTPUT_PRECISE_IL.sh

```bash
#!/bin/bash
SAMP=$1
TE=$2
bwa mem CK2021.05072022.MASKED.fasta ${SAMP}_${TE}.split_reads_L.fas | \
samtools view -q 1 -F 16 | \
awk '{match($6, /([0-9]+)M/, a); print $3"\t"$4+a[1]-1"\t"$4+a[1]}' > ${SAMP}_${TE}.PRECISE.SPLIT_IL.bed

bwa mem CK2021.05072022.MASKED.fasta ${SAMP}_${TE}.split_reads_L.fas | \
samtools view -q 1 -f 16 | \
awk '{match($6, /([0-9]+)M/, a); print $3"\t"$4-1"\t"$4}' >> ${SAMP}_${TE}.PRECISE.SPLIT_IL.bed

bwa mem CK2021.05072022.MASKED.fasta ${SAMP}_${TE}.split_reads_R.fas | \
samtools view -q 1 -f 16 | \
awk '{match($6, /([0-9]+)M/, a); print $3"\t"$4+a[1]-2"\t"$4+a[1]-1}' >> ${SAMP}_${TE}.PRECISE.SPLIT_IL.bed

bwa mem CK2021.05072022.MASKED.fasta ${SAMP}_${TE}.split_reads_R.fas | \
samtools view -q 1 -F 16 | \
awk '{match($6, /([0-9]+)M/, a); print $3"\t"$4-1"\t"$4}' >> ${SAMP}_${TE}.PRECISE.SPLIT_IL.bed

sort -u ${SAMP}_${TE}.PRECISE.SPLIT_IL.bed | bedtools sort -i - | \
bedtools merge -d 20 -i - > ${SAMP}_${TE}.PRECISE.SPLIT_IL_MERGED.bed

bedtools intersect -c -a ${SAMP}_${TE}.PRECISE.SPLIT_IL_MERGED.bed \
-b ${SAMP}_${TE}.PRECISE.SPLIT_IL.bed > ${SAMP}_${TE}.PRECISE_IL.PRE_SPLIT_COUNT.tsv
bedtools intersect -c -a ${SAMP}_${TE}.PRECISE.SPLIT_IL_MERGED.bed \
-b ${TE}/${SAMP}.${TE}.DISC_SUP.bed > ${SAMP}_${TE}.PRECISE_IL.DISC_COUNT.tsv
bedtools intersect -c -a ${SAMP}_${TE}.PRECISE.SPLIT_IL_MERGED.bed \
-b ${TE}/${SAMP}.${TE}.SPLIT_SUP.bed > ${SAMP}_${TE}.PRECISE_IL.OLD_SPLIT_COUNT.tsv
```

```bash
for TE in TE{1..32} ; do for SAMP in $( cat ALL0510.txt ); do ./OUTPUT_PRECISE_IL.sh $SAMP ${TE} ; done; done
```

### # Filter PRECISE ILMERGE

```bash
for TE in TE{1..32}; do
for SAMP in $( cat ALL0510.txt ); do
awk '$4>0 {print $1"\t"$2"\t"$3}' ${SAMP}_${TE}.PRECISE_IL.DISC_COUNT.tsv | \
bedtools sort -i - | bedtools merge -d 100 -i - > ${SAMP}_${TE}.PRECISE.SPLIT_IL_MERGED_FILT.bed
done
done

for TE in TE{1..32}; do
cat *_${TE}.PRECISE.SPLIT_IL_MERGED_FILT.bed | \
bedtools sort -i - | bedtools merge -i - > ${TE}_ALL.PRECISE.SPLIT_IL_MERGED_FILT.TMP.bed
awk '$3-$2>5' ${TE}_ALL.PRECISE.SPLIT_IL_MERGED_FILT.TMP.bed > ${TE}_ALL.0529_PRECISE_IL.TMP.bed
awk '$3-$2<=5' ${TE}_ALL.PRECISE.SPLIT_IL_MERGED_FILT.TMP.bed | \
bedtools sort -i - | bedtools merge -d 20 -i - >> ${TE}_ALL.0529_PRECISE_IL.TMP.bed
bedtools sort -i ${TE}_ALL.0529_PRECISE_IL.TMP.bed > ${TE}_ALL.0529_PRECISE_IL.TMP1.bed
done

```

### # Correspond Precise ILs to the DVS ILs

```bash
for N in {1..32}; do

# Intact DVS TE member : chr1A:start-end
#bedtools sort -i NEW${N}_INTACT_COPIES.EXT10.bed | bedtools merge -i - | \
#awk '{print $1"\t"$2"\t"$3"\t"$1":"$2+10"-"$3-10}' > TE${N}_INTACT.EXT10.bed

# Unintact DVS TE member : chr1A|start-end
#bedtools intersect -v -a NEW${N}_UNINTACT_COPIES.bed -b TE${N}_INTACT.EXT10.bed **| \
#**awk '{print $1"\t"$2-10"\t"$3+10"\t"$1"|"$2"-"$3}' > TE${N}_UNINTACT_EXT10.bed

#cat TE${N}_INTACT.EXT10.bed TE${N}_UNINTACT_EXT10.bed | \
#bedtools sort -i - > TE${N}_DVS.INTACT_UNINTACT.bed 

uniq TE${N}_DVS.INTACT_UNINTACT.bed > TE${N}_DVS.INTACT_UNINTACT.UNIQ.bed
mv TE${N}_DVS.INTACT_UNINTACT.UNIQ.bed TE${N}_DVS.INTACT_UNINTACT.bed

bedtools intersect -wao -a TE${N}_ALL.0529_PRECISE_IL.TMP1.bed -b TE${N}_DVS.INTACT_UNINTACT.bed | \
awk '$7=="\." {print $1"\t"$2"\t"$3"\tNOT_DVS"} $7!="\." {print $1"\t"$2"\t"$3"\t"$7}' \
> TE${N}.PRECISE_IL.DVS.bed

awk '{print $1":"$2"-"$3"\t"$4}' TE${N}.PRECISE_IL.DVS.bed > TE${N}.PRECISE_IL.DVS_MERGED.tsv
done
```

### # Calculate precise and disc read counts

```bash
for TE in TE{1..32}; do
for SAMP in $( cat ALL0510.txt ); do
bedtools intersect -c -a ${TE}.PRECISE_IL.DVS.bed \
-b ${SAMP}_${TE}.PRECISE.SPLIT_IL.bed > ${SAMP}_${TE}.0529_PRECISE_IL.PRE_SPLIT_COUNT.tsv
bedtools intersect -c -a ${TE}.PRECISE_IL.DVS.bed \
-b ${TE}/${SAMP}.${TE}.DISC_SUP.bed > ${SAMP}_${TE}.0529_PRECISE_IL.DISC_COUNT.tsv
done
done
```

### # Merge all sample read counts

```bash
for TE in TE{1..32}; do
awk '$4~/:/ {print $1":"$2"-"$3"\tINTACT"} 
$4~/\|/ {print $1":"$2"-"$3"\tNONINTACT"} 
$4~/DVS/ {print $1":"$2"-"$3"\t"$4}' ${TE}.PRECISE_IL.DVS.bed > ${TE}.TMP.tsv

for SAMP in $( cat ALL0510.txt ); do
cut -f 5 ${SAMP}_${TE}.0529_PRECISE_IL.PRE_SPLIT_COUNT.tsv | paste ${TE}.TMP.tsv - > ${TE}.TMP1.tsv
mv ${TE}.TMP1.tsv ${TE}.TMP.tsv
done
cat header.txt ${TE}.TMP.tsv > ${TE}_ALL.0529_PRECISE_IL.PRE_SPLIT_COUNT.tsv

awk '$4~/:/ {print $1":"$2"-"$3"\tINTACT"} 
$4~/\|/ {print $1":"$2"-"$3"\tNONINTACT"} 
$4~/DVS/ {print $1":"$2"-"$3"\t"$4}' ${TE}.PRECISE_IL.DVS.bed > ${TE}.TMP.tsv

awk 'BEGIN {print "name\tDVS"} {print $1":"$2"-"$3"\t"$4}' ${TE}.PRECISE_IL.DVS.bed > ${TE}.PRECISE_IL.DVS.tsv

for SAMP in $( cat ALL0510.txt ); do
cut -f 5 ${SAMP}_${TE}.0529_PRECISE_IL.DISC_COUNT.tsv | paste ${TE}.TMP.tsv - > ${TE}.TMP1.tsv
mv ${TE}.TMP1.tsv ${TE}.TMP.tsv
done
cat header.txt ${TE}.TMP.tsv > ${TE}_ALL.0529_PRECISE_IL.DISC_COUNT.tsv
done
```

```bash
TE=TE32
for TE in TE{1..15}; do
cut -f 2 ${TE}_GENOTYPE_FILT.tsv | sed 1d | sed "s/:/\t/g;s/-/\t/g" > ${TE}_GENOTYPE_FILT.bed
cut -f 2 ${TE}_GENOTYPE_FILT.tsv | sed 1d > TEMP.tsv
for SAMP in $( cat ALL0510.txt ); do
bedtools intersect -c -a ${TE}_GENOTYPE_FILT.bed -b ${SAMP}_${TE}.PRECISE.SPLIT_IL.bed | cut -f 4 | paste TEMP.tsv - > TEMP1.tsv
mv TEMP1.tsv TEMP.tsv
done
cat header.txt TEMP.tsv > ${TE}_GENOTYPE_FILT.PRECISE_SPLIT_READ_COUNT.tsv
cut -f 2 ${TE}_GENOTYPE_FILT.tsv | sed 1d > TEMP.tsv
for SAMP in $( cat ALL0510.txt ); do
bedtools intersect -c -a ${TE}_GENOTYPE_FILT.bed -b ${TE}/${SAMP}.${TE}.DISC_SUP.bed | cut -f 4 | paste TEMP.tsv - > TEMP1.tsv
mv TEMP1.tsv TEMP.tsv
done
cat header.txt TEMP.tsv > ${TE}_GENOTYPE_FILT.DISC_READ_COUNT.tsv
done

```

### # NEW_SUM_SPLIT_DISC_COUNT.py

```python
#!/usr/bin/env python3
import pandas as pd
import sys

TE=sys.argv[1]
SPLIT=pd.read_table(f'{TE}_ALL.0529_PRECISE_IL.PRE_SPLIT_COUNT.tsv', header=0)
SPLIT.drop(['DVS','Unnamed: 231'], axis=1, inplace=True)
DISC=pd.read_table(f'{TE}_ALL.0529_PRECISE_IL.DISC_COUNT.tsv', header=0)
DISC.drop(['DVS','Unnamed: 231'], axis=1, inplace=True)
IL_DVS=pd.read_table(f'{TE}.PRECISE_IL.DVS.tsv', header=0)
SPLIT=SPLIT.merge(IL_DVS,on='name')
DISC=DISC.merge(IL_DVS,on='name')
SPLIT_DVS=SPLIT.loc[~SPLIT['DVS'].str.contains('DVS')]
SPLIT_NONDVS=SPLIT.loc[SPLIT['DVS'].str.contains('DVS')]
DISC_DVS=DISC.loc[~DISC['DVS'].str.contains('DVS')]
DISC_NONDVS=DISC.loc[DISC['DVS'].str.contains('DVS')]
SPLIT_DVS = SPLIT_DVS.set_index('name').groupby('DVS').sum().reset_index()
DISC_DVS = DISC_DVS.set_index('name').groupby('DVS').sum().reset_index()
SPLIT_DVS['name']=SPLIT_DVS['DVS']
DISC_DVS['name']=DISC_DVS['DVS']
SPLIT=pd.concat([SPLIT_DVS,SPLIT_NONDVS])
DISC=pd.concat([DISC_DVS,DISC_NONDVS])

SPLIT.to_csv(f'ALL_SAMP.{TE}_GENOTYPE_0529.SPLIT_READ_COUNT.Unfilt.tsv', sep="\t", header=True, index=False)
DISC.to_csv(f'ALL_SAMP.{TE}_GENOTYPE_0529.DISC_READ_COUNT.Unfilt.tsv', sep="\t", header=True, index=False)

SPLIT.set_index(['name', 'DVS'], inplace=True)
DISC.set_index(['name', 'DVS'], inplace=True)

SUM=SPLIT+DISC
SUM.reset_index().to_csv(f'ALL_SAMP.{TE}_GENOTYPE_0529.SPLIT_DISC_READ_SUM.Unfilt.tsv', sep="\t", header=True, index=False)

SPLIT_MASK = SPLIT.mask(SPLIT>0, 1)
DISC_MASK = DISC * SPLIT.mask(SPLIT>0, 1)
SPLIT_MASK = SPLIT * DISC.mask(DISC>0, 1)

SUM=SPLIT_MASK+DISC_MASK
SUM=SUM.mask(SUM < 2, 0)
SUM=SUM.mask(SUM > 1, 1)
SUM.reset_index(inplace=True)
SUM.to_csv(f'{TE}_0529_GENOTYPE.tsv', sep="\t", header=True, index=False)

SUM=SPLIT_MASK+DISC_MASK
SUM=SUM.mask(SUM < 2, 0)
SUM.reset_index(inplace=True)
SUM.to_csv(f'ALL_SAMP.{TE}_GENOTYPE_0529.SPLIT_DISC.READ_COUNT_SUM.tsv', sep="\t", header=True, index=False)
```

### # SUM_SPLIT_DISC_COUNT_06062022.py

```python
#!/usr/bin/env python3
import pandas as pd
import sys

TE=sys.argv[1]
SPLIT=pd.read_table(f'{TE}_ALL.0529_PRECISE_IL.PRE_SPLIT_COUNT.tsv', header=0)
SPLIT.drop(['DVS','Unnamed: 231'], axis=1, inplace=True)
DISC=pd.read_table(f'{TE}_ALL.0529_PRECISE_IL.DISC_COUNT.tsv', header=0)
DISC.drop(['DVS','Unnamed: 231'], axis=1, inplace=True)
IL_DVS=pd.read_table(f'{TE}.PRECISE_IL.DVS.tsv', header=0)
SPLIT=SPLIT.merge(IL_DVS,on='name')
DISC=DISC.merge(IL_DVS,on='name')
SPLIT_DVS=SPLIT.loc[~SPLIT['DVS'].str.contains('DVS')]
SPLIT_NONDVS=SPLIT.loc[SPLIT['DVS'].str.contains('DVS')]
DISC_DVS=DISC.loc[~DISC['DVS'].str.contains('DVS')]
DISC_NONDVS=DISC.loc[DISC['DVS'].str.contains('DVS')]
SPLIT_DVS = SPLIT_DVS.set_index('name').groupby('DVS').sum().reset_index()
DISC_DVS = DISC_DVS.set_index('name').groupby('DVS').sum().reset_index()
SPLIT_DVS['name']=SPLIT_DVS['DVS']
DISC_DVS['name']=DISC_DVS['DVS']
SPLIT=pd.concat([SPLIT_DVS,SPLIT_NONDVS])
DISC=pd.concat([DISC_DVS,DISC_NONDVS])

SPLIT.to_csv(f'ALL_SAMP.{TE}_GENOTYPE_0529.SPLIT_READ_COUNT.Unfilt.tsv', sep="\t", header=True, index=False)
DISC.to_csv(f'ALL_SAMP.{TE}_GENOTYPE_0529.DISC_READ_COUNT.Unfilt.tsv', sep="\t", header=True, index=False)

SPLIT.set_index(['name', 'DVS'], inplace=True)
DISC.set_index(['name', 'DVS'], inplace=True)

SUM=SPLIT+DISC
SUM.reset_index().to_csv(f'ALL_SAMP.{TE}_GENOTYPE_0529.SPLIT_DISC_READ_SUM.Unfilt.tsv', sep="\t", header=True, index=False)

SUM = SUM.mask(SUM<2, 0)
SUM.reset_index().to_csv(f'ALL_SAMP.{TE}_GENOTYPE_0606.SPLIT_DISC.FILT_READ_COUNT_SUM.tsv', sep="\t", header=True, index=False)
#SPLIT_MASK = SPLIT.mask(SPLIT>0, 1)
#DISC_MASK = DISC * SPLIT.mask(SPLIT>0, 1)
#SPLIT_MASK = SPLIT * DISC.mask(DISC>0, 1)

#SUM=SPLIT_MASK+DISC_MASK
#SUM=SUM.mask(SUM < 2, 0)
SUM=SUM.mask(SUM > 1, 1)
#SUM.reset_index(inplace=True)
SUM.reset_index().to_csv(f'{TE}_0606_GENOTYPE.tsv', sep="\t", header=True, index=False)

#SUM=SPLIT_MASK+DISC_MASK
#SUM=SUM.mask(SUM < 2, 0)
#SUM.reset_index(inplace=True)
#SUM.to_csv(f'ALL_SAMP.{TE}_GENOTYPE_0529.SPLIT_DISC.READ_COUNT_SUM.tsv', sep="\t", header=True, index=False)
```

### # Merge old and REST MASK_W_count and W_count

```bash
# Correspond new ILs with the old ILs
for TE in TE{1..32}; do 
bedtools intersect -wo -a GENOTYPE/${TE}_0529_GENOTYPE_NONDVS.bed -b GENOTYPE/${TE}_GENOTYPE.bed | \
awk '!a[$4]++ {print $4"\t"$8}' > ${TE}.GT0529_OLDIL.tsv 
done

for TE in TE{1..32} ; do cat ALL0510.txt | parallel -j 20 ./MERGE_OLD_REST_W_COUNT.py {} $TE ; done

for TE in TE{1..32} ; do
rm -f ${TE}.0529IL_ALLELE_MASKED.txt
cut -f 4 ${TE}.A_IL_BMASKED.tsv | sed "s/:/\t/g;s/-/\t/g" | awk '{print $1"\t"$2+500"\t"$3-500}' | bedtools intersect -wo -a GENOTYPE/${TE}_0529_GENOTYPE_NONDVS.bed -b - | cut -f 4 >> ${TE}.0529IL_ALLELE_MASKED.txt
cut -f 4 ${TE}.B_IL_AMASKED.tsv | sed "s/:/\t/g;s/-/\t/g" | awk '{print $1"\t"$2+500"\t"$3-500}' | bedtools intersect -wo -a GENOTYPE/${TE}_0529_GENOTYPE_NONDVS.bed -b - | cut -f 4 >> ${TE}.0529IL_ALLELE_MASKED.txt
cut -f 4 ${TE}_REST.A_IL_BMASKED.tsv | sed "s/:/\t/g;s/-/\t/g" | awk '{print $1"\t"$2+500"\t"$3-500}' | bedtools intersect -wo -a GENOTYPE/${TE}_0529_GENOTYPE_NONDVS.bed -b - | cut -f 4 >> ${TE}.0529IL_ALLELE_MASKED.txt
cut -f 4 ${TE}_REST.B_IL_AMASKED.tsv | sed "s/:/\t/g;s/-/\t/g" | awk '{print $1"\t"$2+500"\t"$3-500}' | bedtools intersect -wo -a GENOTYPE/${TE}_0529_GENOTYPE_NONDVS.bed -b - | cut -f 4 >> ${TE}.0529IL_ALLELE_MASKED.txt
sort -u ${TE}.0529IL_ALLELE_MASKED.txt > ${TE}.0529IL_ALLELE_MASKED.txt1
mv ${TE}.0529IL_ALLELE_MASKED.txt1 ${TE}.0529IL_ALLELE_MASKED.txt
done

```

 ###  Content of MERGE_OLD_REST_W_COUNT.py

```python
#!/usr/bin/env python3

import pandas as pd
import sys

SAMP=sys.argv[1]
TE=sys.argv[2]

for COUNT in ['MASK_W', 'W']:
    old=pd.read_table(f'W_COUNT/{SAMP}_{TE}.{COUNT}_count.tsv',names=['old',f'{COUNT}_COUNT'])
    GT0529_OLD = pd.read_table(f'{TE}.GT0529_OLDIL.tsv', names=['name','old'])
    REST = pd.read_table(f"REST_MASK_W_COUNT/{SAMP}_{TE}_REST.{COUNT}_count.tsv",names=['name', f'{COUNT}_COUNT'])
    OLD=old.merge(GT0529_OLD,on='old')[['name',f'{COUNT}_COUNT']]
    pd.concat([OLD,REST],ignore_index=True).sort_values(by=['name']).to_csv(f'{SAMP}_{TE}.OLD_REST_MERGED.{COUNT}_count.tsv',header=True,index=False,sep="\t")
```

### # NEW_Calculate_MRATIO.py

```python
#!/usr/bin/env python3

import pandas as pd
import sys

TE=sys.argv[1]

MASK = pd.read_table(f'{TE}.0529IL_ALLELE_MASKED.txt',header=None)[0].drop_duplicates()

M_COUNT = pd.read_table(f'ALL_SAMP.{TE}_GENOTYPE_0529.SPLIT_DISC.READ_COUNT_SUM.tsv', header=0)
#FILT=M_COUNT['name'].tolist()
#M_COUNT = M_COUNT.loc[M_COUNT['name'].isin(WMERGE_COUNT['name'])]
M_COUNT_MASK = M_COUNT.loc[M_COUNT['name'].isin(MASK)].copy()
M_COUNT_NOTMASK = M_COUNT.loc[~M_COUNT['name'].isin(MASK)].copy()
M_COUNT_MASK.set_index(['name','DVS'], inplace=True)
M_COUNT_NOTMASK.set_index(['name','DVS'], inplace=True)

W_COUNT = pd.read_table(f'GENOTYPE/{TE}_0529_GENOTYPE.tsv', header=0)
MASKW_COUNT = pd.read_table(f'GENOTYPE/{TE}_0529_GENOTYPE.tsv', header=0)
SAMP = pd.read_table('ALL0510.txt', header=None)

for SRR in SAMP[0]:
    W_COUNT[SRR] = 0
    W = pd.read_table(f'0529_W_count/{SRR}_{TE}.OLD_REST_MERGED.W_count.tsv', header=0)
    for IL in W['name']:
        if IL in W_COUNT['name'].values:
            W_COUNT.loc[W_COUNT['name']==IL,SRR] = W.loc[W['name']==IL,'W_COUNT'].to_list()[0]

for SRR in SAMP[0]:
    MASKW_COUNT[SRR] = 0
    MASKW = pd.read_table(f'0529_W_count/{SRR}_{TE}.OLD_REST_MERGED.MASK_W_count.tsv', header=0)
    for IL in MASKW['name']:
        if IL in MASKW_COUNT['name'].values:
            MASKW_COUNT.loc[MASKW_COUNT['name']==IL,SRR] = MASKW.loc[MASKW['name']==IL,'MASK_W_COUNT'].to_list()[0]

W_COUNT = W_COUNT.set_index(['name','DVS'])
MASKW_COUNT = MASKW_COUNT.set_index(['name','DVS'])
WMERGE_COUNT = W_COUNT + MASKW_COUNT
WMERGE_COUNT.to_csv(f'{TE}_0529_MERGED_W_READ_COUNT.tsv', sep="\t", header=True, index=True)
WMERGE_COUNT = WMERGE_COUNT.reset_index()
WMERGE_COUNT_MASK = WMERGE_COUNT.loc[WMERGE_COUNT['name'].isin(MASK)].copy()
WMERGE_COUNT_NOTMASK = WMERGE_COUNT.loc[~WMERGE_COUNT['name'].isin(MASK)].copy()
WMERGE_COUNT_MASK.set_index(['name','DVS'], inplace=True)
WMERGE_COUNT_NOTMASK.set_index(['name','DVS'], inplace=True)

NOTMASK_MRATIO = M_COUNT_NOTMASK/2/(M_COUNT_NOTMASK/2+WMERGE_COUNT_NOTMASK+0.001)
NOTMASK_MRATIO = NOTMASK_MRATIO.mask(NOTMASK_MRATIO>1,1)
NOTMASK_MRATIO.reset_index(inplace=True)
MASK_MRATIO = M_COUNT_MASK/(M_COUNT_MASK/2+WMERGE_COUNT_MASK+0.001)
MASK_MRATIO = MASK_MRATIO.mask(MASK_MRATIO>1,1)
MASK_MRATIO.reset_index(inplace=True)

pd.concat([NOTMASK_MRATIO, MASK_MRATIO], ignore_index=True).to_csv(f'{TE}_0529_FILT_MRATIO.tsv', sep="\t", header=True, index=False)
```

[Group IL count analysis](https://www.notion.so/Group-IL-count-analysis-c9f79f8da4ba415881fb4a9df92d3038?pvs=21)

[K-Mer Statistics](https://www.notion.so/K-Mer-Statistics-27c830f91108474dab9a9ffa60d63e71?pvs=21)

### # Filter un-intact DVS TEs and nested TE ILs

```bash
for TE in TE{1..32} ; do 
awk 'NR==1 {print $0} $2!~/\|/ && $1~/chr/ {print $0}' ${TE}_0606_GENOTYPE.tsv > ${TE}_0606_GENOTYPE_UNINTACT_FILT.tsv
done
```

[IL-based phylogeny](https://www.notion.so/IL-based-phylogeny-cd5996a7754742c6afe1d14c868c44ae?pvs=21)

[Scanning ILs in CCS of two Valencia](https://www.notion.so/Scanning-ILs-in-CCS-of-two-Valencia-189dd355e56c47dc87549815002e0d19?pvs=21)

[Transposon in evolution](https://www.notion.so/Transposon-in-evolution-9a18176c0ae64d9bbdd8b2ad4377bf00?pvs=21)
