# Sweet orange TE insertion locus (IL) scanning and novel IL detection

## Pipeline Overview
A pipeline for detecting transposon insertion loci (ILs) in sweet orange and related citrus genomes.

---

## Environment Configuration

### 1. Clone the Repository  
**Do not alter the directory structure** to avoid pipeline failures.  
```bash
git clone https://github.com/TheLuoFengLab/Transposon-insertion-locus-simulation-and-scanning.git
chmod 555 -Rf Transposon-insertion-locus-simulation-and-scanning
mv Transposon-insertion-locus-simulation-and-scanning SWOIL  # Optional: rename directory
```

### 2. Install Dependencies via Conda  
```bash
conda env create -f SWOIL/te_detection.yml
```

### 3. Download Reference Genomes  
Download `DVS_reference.zip` from [Figshare](https://doi.org/10.6084/m9.figshare.28737917.v1), then:  
```bash
unzip DVS_reference.zip
mv DVS.fasta SWOIL/REF/
mv DVS.masked.34TE.fasta SWOIL/REF/
```

---

## Running the Pipeline

### Command-Line Help  
```bash
./IL_SCAN.sh -h
```
Output:
```
Usage: SWOIL/IL_SCAN.sh -r REFERENCE -1 READ1 -2 READ2 -o OUT_PREFIX [-b BED] [-t THREADS] [-c CHUNKS]
Transposon insertion detection pipeline

Mandatory options:
  -r, --ref     Reference genome (default: /home/bwu4/IL/REF/DVS.fasta)
  -u, --unmask  Unmasked reference genome (default: /home/bwu4/IL/REF/DVS.masked.34TE.fasta)
  -1, --read1   Read 1 (FASTQ; gz supported)
  -2, --read2   Read 2 (FASTQ; gz supported)
  -o, --out     Output prefix

Optional options:
  -b, --bed     BED file of known ILs (default: /home/bwu4/IL/known_ILs_bed/SWO_ILs_merged.bed)
  -t, --threads Threads for BWA mapping (default: 1)
  -c, --chunks  BAM chunks for parallel processing (default: 1)
```

### Example Execution  
```bash
SWOIL/IL_SCAN.sh -1 CK_1.clean.fastq.gz -2 CK_2.clean.fastq.gz -o ${PREFIX} -t 12 -c 12
```
**Output file:** `${PREFIX}.SWOTE_FINAL.tsv`

---

## Output Format

### Example Output  
```plaintext
Chrom   Start   End     TEfamily        Type    TypeA   TypeB   TypeC   PASS_FILTER
DVS1A   6449981 6456365 CiMULE1 DVS_specific    20      0       0       Yes
DVS1A   11954575        11960959        CiMULE1 DVS_specific    9       0       0       Yes
...
DVS5A   36494523        36494563        CiGYPSY5        Novel   1       1       0       Yes
```

### Column Descriptions  
| Column         | Description                                                                                  |
|----------------|----------------------------------------------------------------------------------------------|
| **Chrom**      | Chromosome ID                                                                                |
| **Start**      | IL start coordinate (exact IL left - 20 bp)                                                 |
| **End**        | IL end coordinate (exact IL right + 20 bp)                                                  |
| **TEfamily**   | TE family name                                                                               |
| **Type**       | IL type: `DVS_specific`, `modern`, `Valencia`, `Novel`, etc.                                 |
| **TypeA**      | Split-mapped read count (max 1 per read pair)                                               |
| **TypeB**      | Discordant read pair count                                                                   |
| **TypeC**      | Concordant read pairs supporting wild-type                                                   |
| **PASS_FILTER**| Passed if: 1. `TypeA > 0` and `TypeA + TypeB ≥ 2` 2. `TypeA + TypeB + TypeC < 2.5× sequencing depth` |

---

## Notes  
- **Reference Files:** Both `DVS.fasta` and `DVS.masked.34TE.fasta` must reside in `SWOIL/REF/`.
- **Novel ILs:** Marked as `Novel` in the `Type` column.
- **Filtering:** `PASS_FILTER=Yes` indicates ILs passing depth and support thresholds.
