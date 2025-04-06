#!/bin/bash
set -eo pipefail

DIR=$(dirname "$(realpath "$0")")
SCRIPT_DIR=${DIR}/scripts

# Help information
usage() {
    echo "Usage: $0 -r REFERENCE -1 READ1 -2 READ2 -o OUT_PREFIX [-b BED] [-t THREADS] [-c CHUNKS]"
    echo "Transposon insertion detection pipeline"
    echo ""
    echo "Mandatory options:"
    echo "  -r, --ref     Reference genome"
    echo "  -u, --unmask  Unmasked reference genome"
    echo "  -1, --read1   First read file (FASTQ)"
    echo "  -2, --read2   Second read file (FASTQ)"
    echo "  -o, --out     Output prefix"
    echo ""
    echo "Optional options:"
    echo "  -b, --bed     BED file with known ILs (optional)"
    echo "  -t, --threads Number of threads (default: 1)"
    echo "  -c, --chunks  Number of BAM chunks (default: 1)"
    echo "  -h, --help    Show this help message"
    exit 0
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -r|--ref) REFERENCE="$2"; shift 2;;
        -u|--unmask) UNMASK_REF="$2"; shift 2;;
        -1|--read1) READ1="$2"; shift 2;;
        -2|--read2) READ2="$2"; shift 2;;
        -o|--out) OUT_PREFIX="$2"; shift 2;;
        -b|--bed) BED_FILE="$2"; shift 2;;
        -t|--threads) THREADS="$2"; shift 2;;
        -c|--chunks) CHUNKS="$2"; shift 2;;
        -h|--help) usage;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

# Check mandatory arguments
if [[ -z $READ1 || -z $READ2 || -z $OUT_PREFIX ]]; then
    echo "Missing required arguments!"
    usage
    exit 1
fi

# Set defaults
THREADS=${THREADS:-1}
CHUNKS=${CHUNKS:-1}
SAMTOOLS_THREADS=$(( THREADS < 4 ? THREADS : 4 ))
REFERENCE=${REFERENCE:-${DIR}/REF/DVS.masked.34TE.fasta}
UNMASK_REF=${UNMASK_REF:-${DIR}/REF/DVS.fasta}
BED_FILE=${BED_FILE:-${DIR}/known_ILs_bed/SWO_ILs_merged.bed}

if [[ ! -f $REFERENCE ]]; then
    echo "${REFERENCE} missing! Quit running."
    usage
    exit 1
fi

if [[ ! -f $UNMASK_REF || -z $BED_FILE ]]; then
    echo "${UNMASK_REF} missing! Quit running."
    usage
    exit 1
fi

if [[ ! -f $BED_FILE ]]; then
    echo "${BED_FILE} missing! Quit running."
    usage
    exit 1
fi

# Step 1: Alignment
STEP1_DONE="${OUT_PREFIX}.alignment.done"
STEP1_BAM="${OUT_PREFIX}_namesorted.bam"
STEP1_UNMASK_BAM="${OUT_PREFIX}.unmasked_pos_srt.bam"

if [[ (! -f "$STEP1_DONE") || (! -f "${OUT_PREFIX}.split_reads.txt") || (! -f "$STEP1_BAM") ]]; then
    echo "Step 1/4: Aligning and sorting reads"
    rm -f "$STEP1_DONE"

    # Index reference if needed
    if [[ ! -f "${UNMASK_REF}.bwt" ]]; then
        echo "Indexing reference genome..."
        bwa index "$UNMASK_REF" 2>/dev/null
    fi
    if [[ ! -f "$STEP1_UNMASK_BAM" ]]; then
        bwa mem -T 20 -Y -t 12 $UNMASK_REF $READ1 $READ2 | samtools sort -@ 4 -o $STEP1_UNMASK_BAM
    fi

    samtools view -@ $THREADS -F 0x900 -d "SA" DVS.unmasked_pos_srt.bam | cut -f 1 | sort -u > ${OUT_PREFIX}.split_reads.txt
    
    # Index reference if needed
    if [[ ! -f "${REFERENCE}.bwt" ]]; then
        echo "Indexing reference genome..."
        bwa index "$REFERENCE" 2>/dev/null
    fi

    # Alignment and sorting
    bwa mem -T 20 -Y -t "$THREADS" "$REFERENCE" "$READ1" "$READ2" | \
    samtools sort -n -@ "$SAMTOOLS_THREADS" -o "$STEP1_BAM" - && touch "$STEP1_DONE"
else
    echo "Step 1/4: Alignment completed, skipping..."
fi

# Step 2: Fixmate and Sorting
STEP2_DONE="${OUT_PREFIX}.fixmate.done"
FIXMATE_BAM="${OUT_PREFIX}_fixmate.bam"
POS_SORTED_BAM="${OUT_PREFIX}_pos_sorted.bam"
DEDUP_BAM="${OUT_PREFIX}_dedup.bam"
NAME_SORTED_BAM="${OUT_PREFIX}_name_sorted.bam"

if [[ ! -f "$STEP2_DONE" ]]; then
    echo "Step 2/4: Processing duplicates"
    if [[ ! -f "$FIXMATE_BAM" ]]; then
        # Fixmate repair
        samtools fixmate -@ "$SAMTOOLS_THREADS" -m "$STEP1_BAM" "$FIXMATE_BAM"
    fi

    if [[ ! -f "$POS_SORTED_BAM" ]]; then
        # Position-sorting required for rmdup [2,5](@ref)
        samtools sort -@ "$SAMTOOLS_THREADS" -o "$POS_SORTED_BAM" "$FIXMATE_BAM"
    fi
    if [[ ! -f  "$DEDUP_BAM" ]]; then
        # PCR duplicate removal (single-end mode)
        samtools rmdup -s "$POS_SORTED_BAM" "$DEDUP_BAM" 2> "${OUT_PREFIX}.rmdup.log"
    fi
    if [[ ! -f "$NAME_SORTED_BAM" ]]; then
        # Re-name-sort for downstream processing [1,4](@ref)
        samtools sort -n -@ "$SAMTOOLS_THREADS" -o "$NAME_SORTED_BAM" "$DEDUP_BAM"
    fi

    touch "$STEP2_DONE"
else
	echo "Step 2/4: Removing duplicates completed, skipping..."
fi

# Step 3: Chunk Preparation
STEP3_DONE="${OUT_PREFIX}.chunk_prep.done"
CHUNK_DIR="${OUT_PREFIX}_chunks"
READ_NAMES="${OUT_PREFIX}_mapq_names.txt"

if [[ (! -f "$STEP3_DONE") ]]; then
    echo "Step 3/4: Preparing BAM chunks"
    rm -rf "$CHUNK_DIR"
    mkdir -p "$CHUNK_DIR"
    
    # Extract mapped read names with MAPQ > 0
    samtools view -F 4 -q 1 -@ "$SAMTOOLS_THREADS" "$NAME_SORTED_BAM" | \
        awk '{print $1}' | sort -u > "$READ_NAMES"
    
    if [[ $CHUNKS -gt 1 ]]; then
        # Split names into balanced chunks
        split -d -n l/"$CHUNKS" "$READ_NAMES" "${CHUNK_DIR}/chunk_"
        
        # Create chunk-specific BAMs
        for CHUNK_FILE in "${CHUNK_DIR}"/chunk_*; do
            samtools view -@ "$SAMTOOLS_THREADS" -N "$CHUNK_FILE" \
                -b -o "${CHUNK_FILE}.bam" "$NAME_SORTED_BAM"
        done
    else
        # Single chunk symlink
        ln -sf "$NAME_SORTED_BAM" "${CHUNK_DIR}/chunk_00.bam"
    fi
    
    touch "$STEP3_DONE"
else
    echo "Step 3/4: Chunk preparation completed, skipping..."
fi

# Step 4: Parallel Detection and Merging
STEP4_DONE="${OUT_PREFIX}.detection.done"
STEP4_PARTS="${OUT_PREFIX}_il_parts"
STEP4_REPORT="${OUT_PREFIX}_il_report_raw.tsv"
MERGED_BED="${OUT_PREFIX}_temp_ils.bed"

if [[ (! -f "$STEP4_DONE") ]]; then
    echo "Step 4/4: Detecting insertion loci"
    rm -rf "$STEP4_PARTS"
    mkdir -p "$STEP4_PARTS"

    # Phase 1: First pass IL detection
    samtools view -@8 -N ${OUT_PREFIX}.split_reads.txt -o ${OUT_PREFIX}_dedup.split.bam ${OUT_PREFIX}_dedup.bam
    python3 ${SCRIPT_DIR}/first_pass.py -b ${OUT_PREFIX}_dedup.split.bam ${BED_FILE:+--bed $BED_FILE} -o ${MERGED_BED}
    
    if [[ $CHUNKS -gt 1 ]]; then
        # Parallel execution
        # Phase 2: Second pass counting
        mkdir -p "${OUT_PREFIX}_secondpass"
        find "$CHUNK_DIR" -name "chunk_*.bam" -print0 | \
            parallel -0 -j "$THREADS" \
                "python3 ${SCRIPT_DIR}/second_pass.py -b {} -m $MERGED_BED -o ${OUT_PREFIX}_secondpass/secondpass_{#}.tsv"

        # Merge final reports
        python3 ${SCRIPT_DIR}/merge_reports.py \
            -i "${OUT_PREFIX}_secondpass"/secondpass_*.tsv \
            -o "$STEP4_REPORT"

        touch "$STEP4_DONE"
    else
        # Single chunk processing
        python3 ${SCRIPT_DIR}/second_pass.py \
            -m $MERGED_BED \
            -b "${OUT_PREFIX}_dedup.bam" \
            -o "$STEP4_REPORT" \
            ${BED_FILE:+--bed "$BED_FILE"}
    fi
    
    samtools depth -@ $THREADS $STEP1_UNMASK_BAM > ${OUT_PREFIX}.depth
    DEPTH=$( shuf -n 1000000 ${OUT_PREFIX}.depth | cut -f 3 | sort -n | uniq -c | sort -nk1,1 | tail -n1 | awk '{print $2}') 
    
    awk -v x=$DEPTH 'BEGIN {OFS="\t"} NR==1 {print $0} $5!~/Ancestral/ && $5!~/Unintact/ && $9=="Yes" && ($6+$7+$8)<2.5*x' $STEP4_REPORT > ${OUT_PREFIX}.temp.tsv
    
    python3 ${SCRIPT_DIR}/id_conversion.py ${OUT_PREFIX}.temp.tsv IL_CONVERSION.tsv ${OUT_PREFIX}.SWOTE_FINAL.tsv

    touch "$STEP4_DONE"

else
    echo "Step 4/4: Detection completed, skipping..."
fi

echo "Pipeline completed. Results in ${OUT_PREFIX}.SWOTE_FINAL.tsv"
