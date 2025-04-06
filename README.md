# Sweet orange TE insertion locus (IL) scanning and novel IL detection

## Pipeline environment configuration
### Download the repository. Do not alter the dir structure otherwise the running of IL_SCAN.sh will fail due to not finding the target file.
  git clone https://github.com/TheLuoFengLab/Transposon-insertion-locus-simulation-and-scanning.git
  chmod 555 -Rf Transposon-insertion-locus-simulation-and-scanning
  mv Transposon-insertion-locus-simulation-and-scanning SWOIL  # Rename the dir to any name you want

### Install tools in TESCAN.yml
  conda env create -f SWOIL/te_detection.yml

### Download reference files from Figshare
###### #The link for original and modified diploid Valencia sweet orange genomes are compressed in "DVS_reference.zip" available from https://doi.org/10.6084/m9.figshare.28737917.v1.
###### #Both DVS.fasta and DVS.masked.34TE.fasta are required to run the pipeline. And they need to be put in the REF dir.
  unzip DVS_reference.zip
  mv DVS.fasta SWOIL/REF/
  mv DVS.masked.34TE.fasta SWOIL/REF/

### Help information of IL_SCAN.sh
$ ./IL_SCAN.sh -h
Usage: SWOIL/IL_SCAN.sh -r REFERENCE -1 READ1 -2 READ2 -o OUT_PREFIX [-b BED] [-t THREADS] [-c CHUNKS]
Transposon insertion detection pipeline

Mandatory options:
  -r, --ref     Reference genome default as /home/bwu4/IL/REF/DVS.fasta. Usually no need to modify.
  -u, --unmask  Unmasked reference genome default as /home/bwu4/IL/REF/DVS.masked.34TE.fasta. Usually no need to modify.
  -1, --read1   First read file (FASTQ) could either be compressed in gz or not
  -2, --read2   Second read file (FASTQ) could either be compressed in gz or not
  -o, --out     Output prefix

Optional options:
  -b, --bed     BED file with known ILs default as /home/bwu4/IL/known_ILs_bed/SWO_ILs_merged.bed. Usually no need to modify.
  -t, --threads Number of threads to be used for bwa mapping (default: 1)
  -c, --chunks  Number of BAM chunks to be processed parallely. Usually do not exceed core number (default: 1)
  -h, --help    Show this help message

## Running IL_SCAN.sh on pair-end NGS of citrus (working for sweet orange, pummelo, and mandarins)
###### #Input fastq file could either be compressed or uncompressed
SWOIL/IL_SCAN.sh -1 CK_1.clean.fastq.gz -2 CK_2.clean.fastq.gz -o ${PREFIX} -t 12 -c 12  # The final report will be named as ${PREFIX}.SWOTE_FINAL.tsv 

## Examplary output
Chrom   Start   End     TEfamily        Type    TypeA   TypeB   TypeC   PASS_FILTER
DVS1A   6449981 6456365 CiMULE1 DVS_specific    20      0       0       Yes
DVS1A   11954575        11960959        CiMULE1 DVS_specific    9       0       0       Yes
DVS1A   23495527        23501911        CiMULE1 DVS_specific    10      0       0       Yes
DVS6A   269266  269314  CiMULE2 Novel   6       1       12      Yes
DVS8A   13623853        13630816        CiMULE2 non-China       18      0       0       Yes
DVS8A   28575053        28582016        CiMULE2 Valencia        8       0       0       Yes
DVS8B   11997659        12002756        CiMULE2 modern  13      0       0       Yes
DVS9A   3577018 3577066 CiMULE2 Novel   4       1       5       Yes
DVS1B   7382    7422    CiGYPSY1        Novel   4       0       10      Yes
DVS3B   16860726        16863770        CiGYPSY1        modern  9       0       0       Yes
DVS8A   11648570        11648610        CiGYPSY1        Novel   3       0       12      Yes
DVS1B   2608925 2613430 CiHAT1  modern  5       0       0       Yes
DVS3A   18607322        18607362        CiHAT1  Novel   1       1       17      Yes
DVS7A   3975131 3975171 CiHAT1  Novel   1       1       27      Yes
DVS3A   14144040        14144080        CiGYPSY5        Novel   5       2       11      Yes
DVS5A   16481477        16481517        CiGYPSY5        Novel   3       0       13      Yes
DVS5A   36494523        36494563        CiGYPSY5        Novel   1       1       0       Yes
......

Explanations for the columns:
Chrom: chromosome ID
Start: reported chromosome coordinate start for a TE insertion locus (IL), usually reported as exact IL left coordinate - 20 bp
End: reported chromosome coordinate end for a TE insertion locus (IL), usually reported as exact IL right coordinate + 20 bp
TEfamily: name of the TE family contributing to the IL
Type: DVS_specific: ILs specificly detected in the DVS assembly; tag ILs of cultivar groups: modern, non-China, Valencia, Navel, CQ/Jincheng, CQ, et al.; Novel, novel ILs detected in the accession;
TypeA: count of split-mapped reads supporting the IL. A read pair is counted at most once.
TypeB: count of discordantly mapped read pairs supporting the IL. 
TypeC: count of concordantly mapped read pairs across the ILs that support the wild (Reference) type.
PASS_FILTER: ILs that have passed two filters: (1) TypeA>0 and TypeA + TypeB>=2 (2) TypeA + TypeB + TypeC < 2.5-fold sequencing depth (to remove false positives in repeat genome regions) 

