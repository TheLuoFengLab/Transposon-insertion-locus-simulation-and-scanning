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
  awk '$5==1 && $4-$3<=10 && $11-$10>50 && $11-$10<20000 {print $2"\t"$3"\t"$4}' $VAR > ${VAR%.fasta.var.txt}.INS.bed
  # Extracts deletion (DEL) and insertion (INS) variants shorter than 20000 bp and saves them to respective bed files

  cat ${VAR%.fasta.var.txt}.DEL.tsv ${VAR%.fasta.var.txt}.INS.tsv | bedtools sort -i - > ${VAR%.fasta.var.txt}.INDEL.bed
  # Concatenates the deletion and insertion files into a single INDEL.bed file and sort according to the chromosomes and coordinates
done
cat *.INDEL.bed | bedtools merge -i -d

### Step 3: 'for' loop - Intersection analysis with bedtools:
for BED in [list of INDEL.bed files]; do
  bedtools intersect -f 0.95 -F 0.95 -c -a ALL10.INDEL.bed -b $BED > ALL10.${BED}
  # Uses bedtools intersect to find overlaps between the ALL10.INDEL.bed file and each BED file
done

### Step 4: 'for' loop - Creating a combined file:
for BED in ALL10.*.INDEL.bed ; do 
  cut -f 4 $BED | paste TEMP - > TEMP1 ; mv TEMP1 TEMP; 
done
# Iterates over ALL10.*.INDEL.bed files, extracts the 4th field, and merges it into a temporary file
