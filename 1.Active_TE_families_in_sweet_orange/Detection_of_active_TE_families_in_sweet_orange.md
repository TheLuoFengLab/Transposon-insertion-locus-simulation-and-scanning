# Detection of active TE families in sweet orange
1.1 Download sweet orange assemblies
Assemblies information:
|Acession name	| Description	| NCBI accession No. |
| ------------- | ----------- | ------------------ |
| DVS	| Valencia sweet orange from Florida | GCA_022201045.1 (DVS_A), GCA_022201065.1 (DVS_B) |
| T19	| Irradiated Valencia sweet orange mutant T19	| SRR15706502
| T78	| Irradiated Valencia sweet orange mutant T78	| SRR15706505
| HSO	| Di-haploid Valencia sweet orange	| GCA_018105775.1
| VAL	| Valencia sweet orange	| GCA_018104345.1
| NHE	| Newhall navel orange	| GCA_019144195.1
| NW	| Lanlate: a late maturing mutant of navel orange	| GCA_019144185.1
| SO3	| Valencia sweet orange	| GCA_019143665.1
| TCPS1	| a high acid orange from Hunan Province in South China	| GCA_019144155.1
| UKXC	| a local sweet orange variety from Hunan province |	GCA_019144245.1
| BT2	| Bingtangcheng No.2	| GCA_019144225.1

Download assemblies from NCBI:
```bash
# Make sure no useful ncbi_dataset.zip and ncbi_dataset/ are present in current working dir
# The download assemblies are saved in ./assemblies
mkdir assemblies
for ACC in GCA_022201045.1 GCA_022201065.1 GCA_018105775.1 GCA_018104345.1 GCA_019144195.1 GCA_019144185.1 GCA_019143665.1 GCA_019144155.1 GCA_019144245.1 GCA_019144225.1 ; do
  rm -Rf ncbi_dataset ncbi_dataset.zip
    if [[ ! -a assemblies/${ACC}.fasta ]]; then
    datasets download genome accession ${ACC} --include genome
    unzip -o ncbi_dataset.zip
    mv ncbi_dataset/data/${ACC}/${ACC}*.fna assemblies/${ACC}.fasta
    fi
done
```
Download T19 and T78 diploid assemblies
'''
wget --no-check-certificate -O T78.asem.fasta.gz  https://figshare.com/ndownloader/files/45090274 && gzip -d T78.asem.fasta.gz
mv T78.asem.fasta assemblies/
wget --no-check-certificate -O T19.asem.fasta.gz  https://figshare.com/ndownloader/files/45090271 && gzip -d T19.asem.fasta.gz
mv T19.asem.fasta assemblies/
'''

Concatenate DVS_A and DVS_B as into diploid DVS assembly
```
cat assemblies/GCA_022201045.1.fasta assemblies/GCA_022201065.1.fasta > assemblies/DVS.fasta
```


1.2 Align 10 SWO assemblies to DVS and call large indels
```bash
cd assemblies
for FAS in GCA_019144245.1.fasta GCA_019144225.1.fasta GCA_019144195.1.fasta GCA_019144185.1.fasta GCA_019144155.1.fasta \
GCA_019143665.1.fasta GCA_018104345.1.fasta GCA_018105775.1.fasta T78.asem.fasta SF.asem.fasta ; do
  rm -f asm.paf asm.srt.paf
  minimap2 -cx asm5 -t8 --cs DVS.fasta $FAS > asm.paf  
  sort -k6,6 -k8,8n asm.paf > asm.srt.paf             # sort by reference start coordinate
  k8 paftools.js call asm.srt.paf > ${FAS}.var.txt
  awk '$5==1 && $4-$3>50 && $4-$3<20000 && $11-$10<=10 {print $2"\t"$3"\t"$4}' ${FAS}.var.txt > ${FAS%.fasta}.DEL.tsv
  awk '$5==1 && $4-$3<=10 && $11-$10>50 && $11-$10<20000 {print $2"\t"$3"\t"$4}' ${FAS}.var.txt > ${FAS%.fasta}.INS.tsv
  cat ${FAS%.fasta}.DEL.tsv ${FAS%.fasta}.INS.tsv > ${FAS%.fasta}.INDEL.bed
done
```

1.3 

for BED in *.INDEL.bed; do
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
