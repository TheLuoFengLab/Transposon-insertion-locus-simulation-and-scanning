### Step 1: 'for' loop - Processing FASTA files with minimap2, sorting and calling variants:

for FAS in GCA_019144245.1.fasta GCA_019144225.1.fasta GCA_019144195.1.fasta GCA_019144185.1.fasta GCA_019144155.1.fasta GCA_019143665.1.fasta GCA_018104345.1.fasta Csiv4.fasta T78.asem.fasta SF.asem.fasta; do
# Iterates over a list of SWO genome assembly files using DVS diploid genome assembly (CK2021.60.corrected.fasta) as reference

  minimap2 -cx asm5 -t8 --cs CK2021.60.corrected.fasta $FAS > asm.paf  
  # Aligns each FASTA file against a reference genome using minimap2 with specific options (-cx asm5 -t8 --cs)

  sort -k6,6 -k8,8n asm.paf > asm.srt.paf             
  # Sorts the PAF output file by reference start coordinate

  k8 paftools.js call asm.srt.paf > ${FAS}.var.txt
  # Calls variants using paftools.js and saves the output to a .var.txt file for each FASTA file
done

### Step 2:  'for' loop - Extracting and merging deletion and insertion variants:
for VAR in *.fasta.var.txt ; do
# Iterates over all .fasta.var.txt files from previous step
  awk '$5==1 && $4-$3>50 && $4-$3<20000 && $11-$10<=10 {print $2"\t"$3"\t"$4}' $VAR > ${VAR%.fasta.var.txt}.DEL.bed
  awk '$5==1 && $4-$3<=20 && $11-$10>50 && $11-$10<20000 {print $2"\t"$3"\t"$4}' $VAR > ${VAR%.fasta.var.txt}.INS.bed
  # Extracts deletion (DEL) and insertion (INS) variants shorter than 20000 bp and saves them to respective bed files
done

cat *.DEL.bed | bedtools sort -i - > merged.DEL.temp.bed
# Concatenates all .DEL.bed files and sorts them, preparing them for further processing
bedtools intersect -wo -f 0.95 -r -a merged.DEL.temp.bed -b merged.DEL.temp.bed | awk '$2!=$5 || $3!=$6' > merged.temp.DEL.overlaps.bed
# Find all overlaps with over 95% reciprocal overlap within merged.temp.DEL.bed. The awk part filters out self-overlaps (where an interval overlaps with itself).
cut -f 1,2,3 merged.temp.DEL.overlaps.bed | uniq > DEL.duplicate.bed
# Extract unique intervals that are involved in overlaps and stores them in duplicate.DEL.bed.
sort merged.DEL.temp.bed > merged.DEL.temp.srt.bed
sort DEL.duplicate.bed > DEL.duplicate.srt.bed
comm -23 merged.DEL.temp.srt.bed DEL.duplicate.srt.bed > DEL.non-duplicate.bed
# Separate non-duplicate intervals from the original file
while [[ $( wc -l merged.temp.DEL.overlaps.bed | cut -d " " -f 1) -gt 0 ]]; do
awk '$2<=$5 {if ($3>=$6) {print $1"\t"$2"\t"$3} else {print $1"\t"$2"\t"$6}}' merged.temp.DEL.overlaps.bed > temp_duplicate.bed
bedtools intersect -wo -f 0.95 -r -a temp_duplicate.bed -b temp_duplicate.bed | awk '$2!=$5 || $3!=$6' > merged.temp.DEL.overlaps.bed
done
# iteratively merge overlapping regions. This loop continues as long as there are overlaps found in merged.temp.DEL.overlaps.bed. It processes these overlaps and 
# writes them to temp_duplicate.bed, then recalculates overlaps in each iteration.
cat non-duplicate.DEL.bed temp_duplicate.bed | bedtools sort -i - > ALL10.DEL.filtered_overlaps.bed
# Concatenate DEL from all assemblies and remove duplications with reciprocal 95% overlap

cat *.INS.bed | bedtools sort -i - | bedtools merge -d 10 -i - > ALL10.INS.merged.bed
# Concatenate INSs from all assemblies and deduplicate those within 10 bp distances

### Step 3: 'for' loop - Intersection analysis with bedtools:
for BED in Csiv4.DEL.bed T78.asem.DEL.bed SF.asem.DEL.bed GCA_019144245.1.DEL.bed GCA_019144225.1.DEL.bed GCA_019144195.1.DEL.bed GCA_019144185.1.DEL.bed GCA_019144155.1.DEL.bed GCA_019143665.1.DEL.bed GCA_018104345.1.DEL.bed; do
  bedtools intersect -f 0.9 -F 0.9 -c -a ALL10.DEL.filtered_overlaps.bed -b $BED > TEMP.${BED}
  # Uses bedtools intersect to find overlaps between the ALL10.INDEL.bed file and each BED file
done

for BED in *.INS.bed ; do
  bedtools intersect -F 0.99 -c -a ALL10.DEL.filtered_overlaps.bed -b $BED > TEMP.${BED}
done
  

### Step 4: Creating combined files:
echo -e "CHR\tSTART\tEND\t"$( ls TEMP.*.DEL.bed ) | sed "s/ /\t/g" > header.txt
for BED in TEMP.*.DEL.bed ; do 
  cut -f 4 $BED | paste TEMP - > TEMP1
  mv TEMP1 TEMP; 
done
cp header.txt ALL10.DEL.tsv
paste ALL10.DEL.filtered_overlaps.bed TEMP1 >> ALL10.DEL.tsv
# Iterates over ALL10.*.INDEL.bed files, extracts the 4th field, and merges it into a temporary file


