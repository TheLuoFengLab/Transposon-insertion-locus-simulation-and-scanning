#!/bin/bash
set -euo pipefail

DIR=$(dirname "$(realpath "$0")")
SCRIPT_DIR="${DIR}/scripts"

usage() {
    cat <<EOF
Usage: $0 -1 READ1 -2 READ2 -o OUT_PREFIX [options]

Transposon insertion locus detection pipeline.

Mandatory options:
  -1, --read1      First read file, FASTQ/FASTQ.gz
  -2, --read2      Second read file, FASTQ/FASTQ.gz
  -o, --out        Output prefix

Optional options:
  -r, --ref        Masked / TE-augmented scanning reference
                   Default: ${DIR}/REF/DVS.masked.34TE.fasta

  -u, --unmask     Unmasked reference used for split-read discovery and depth
                   Default: ${DIR}/REF/DVS.fasta

  -b, --bed        BED file with known insertion loci
                   Default: ${DIR}/known_ILs_bed/SWO_ILs_merged.bed

  -t, --threads    Number of threads
                   Default: 1

  -c, --chunks     Number of query-name chunks for second_pass.py
                   Default: 1

  -h, --help       Show this help message

Important:
  Chunk BAMs are kept query-name sorted because second_pass.py groups records
  by adjacent read names. Do not coordinate-sort chunk BAMs before second_pass.py.
EOF
    exit 0
}

# -----------------------------
# Parse command-line arguments
# -----------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--ref)
            REFERENCE="$2"
            shift 2
            ;;
        -u|--unmask)
            UNMASK_REF="$2"
            shift 2
            ;;
        -1|--read1)
            READ1="$2"
            shift 2
            ;;
        -2|--read2)
            READ2="$2"
            shift 2
            ;;
        -o|--out)
            OUT_PREFIX="$2"
            shift 2
            ;;
        -b|--bed)
            BED_FILE="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -c|--chunks)
            CHUNKS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            ;;
    esac
done

# -----------------------------
# Required inputs and defaults
# -----------------------------
if [[ -z "${READ1:-}" || -z "${READ2:-}" || -z "${OUT_PREFIX:-}" ]]; then
    echo "ERROR: Missing required arguments." >&2
    usage
fi

THREADS=${THREADS:-1}
CHUNKS=${CHUNKS:-1}
SAMTOOLS_THREADS=$(( THREADS < 4 ? THREADS : 4 ))

REFERENCE=${REFERENCE:-${DIR}/REF/DVS.masked.34TE.fasta}
UNMASK_REF=${UNMASK_REF:-${DIR}/REF/DVS.fasta}
BED_FILE=${BED_FILE:-${DIR}/known_ILs_bed/SWO_ILs_merged.bed}

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    echo "ERROR: --threads must be a positive integer." >&2
    exit 1
fi

if ! [[ "$CHUNKS" =~ ^[0-9]+$ ]] || [[ "$CHUNKS" -lt 1 ]]; then
    echo "ERROR: --chunks must be a positive integer." >&2
    exit 1
fi

# -----------------------------
# Dependency checks
# -----------------------------
for cmd in bwa samtools awk sort split shuf python3; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: Required command not found: $cmd" >&2
        exit 1
    }
done

if [[ "$CHUNKS" -gt 1 ]]; then
    command -v parallel >/dev/null 2>&1 || {
        echo "ERROR: GNU parallel is required when --chunks > 1." >&2
        exit 1
    }
fi

# -----------------------------
# File checks
# -----------------------------
[[ -f "$READ1" ]] || { echo "ERROR: READ1 missing: $READ1" >&2; exit 1; }
[[ -f "$READ2" ]] || { echo "ERROR: READ2 missing: $READ2" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "ERROR: Reference missing: $REFERENCE" >&2; exit 1; }
[[ -f "$UNMASK_REF" ]] || { echo "ERROR: Unmasked reference missing: $UNMASK_REF" >&2; exit 1; }
[[ -f "$BED_FILE" ]] || { echo "ERROR: BED file missing: $BED_FILE" >&2; exit 1; }

[[ -f "${SCRIPT_DIR}/first_pass.py" ]] || { echo "ERROR: Missing ${SCRIPT_DIR}/first_pass.py" >&2; exit 1; }
[[ -f "${SCRIPT_DIR}/second_pass.py" ]] || { echo "ERROR: Missing ${SCRIPT_DIR}/second_pass.py" >&2; exit 1; }
[[ -f "${SCRIPT_DIR}/merge_reports.py" ]] || { echo "ERROR: Missing ${SCRIPT_DIR}/merge_reports.py" >&2; exit 1; }
[[ -f "${SCRIPT_DIR}/id_conversion.py" ]] || { echo "ERROR: Missing ${SCRIPT_DIR}/id_conversion.py" >&2; exit 1; }
[[ -f "${SCRIPT_DIR}/IL_CONVERSION.tsv" ]] || { echo "ERROR: Missing ${SCRIPT_DIR}/IL_CONVERSION.tsv" >&2; exit 1; }

# -----------------------------
# Output file names
# -----------------------------
STEP1_DONE="${OUT_PREFIX}.alignment.done"
STEP1_BAM="${OUT_PREFIX}_namesorted.bam"
STEP1_UNMASK_BAM="${OUT_PREFIX}.unmasked_pos_srt.bam"
SPLIT_READS="${OUT_PREFIX}.split_reads.txt"

STEP2_DONE="${OUT_PREFIX}.fixmate.done"
FIXMATE_BAM="${OUT_PREFIX}_fixmate.bam"
POS_SORTED_BAM="${OUT_PREFIX}_pos_sorted.bam"
DEDUP_BAM="${OUT_PREFIX}_dedup.bam"
NAME_SORTED_BAM="${OUT_PREFIX}_name_sorted.bam"

STEP3_DONE="${OUT_PREFIX}.chunk_prep.done"
STEP3_META="${OUT_PREFIX}.chunk_prep.meta"
CHUNK_DIR="${OUT_PREFIX}_chunks"
READ_NAMES="${OUT_PREFIX}_mapq_names.txt"

STEP4_DONE="${OUT_PREFIX}.detection.done"
STEP4_REPORT="${OUT_PREFIX}_il_report_raw.tsv"
MERGED_BED="${OUT_PREFIX}_temp_ils.bed"
SPLIT_BAM="${OUT_PREFIX}_dedup.split.bam"
SPLIT_POS_BAM="${OUT_PREFIX}_dedup.split.pos.bam"
FINAL_TSV="${OUT_PREFIX}.SWOTE_FINAL.tsv"

# -----------------------------
# Step 1: Alignment
# -----------------------------
if [[ ! -f "$STEP1_DONE" || ! -f "$SPLIT_READS" || ! -f "$STEP1_BAM" ]]; then
    echo "Step 1/4: Alignment"
    rm -f "$STEP1_DONE"

    if [[ ! -f "${UNMASK_REF}.bwt" ]]; then
        echo "Indexing unmasked reference..."
        bwa index "$UNMASK_REF" 2>/dev/null
    fi

    if [[ ! -f "$STEP1_UNMASK_BAM" ]]; then
        echo "Aligning reads to unmasked reference..."
        bwa mem -T 20 -Y -t "$THREADS" "$UNMASK_REF" "$READ1" "$READ2" | \
            samtools sort -@ "$SAMTOOLS_THREADS" -o "$STEP1_UNMASK_BAM" -
    fi

    echo "Extracting SA-tagged split-read names..."
    samtools view -@ "$SAMTOOLS_THREADS" -F 0x900 "$STEP1_UNMASK_BAM" | \
        awk '/SA:Z:/ {print $1}' | \
        sort -u > "$SPLIT_READS"

    if [[ ! -s "$SPLIT_READS" ]]; then
        echo "WARNING: No SA-tagged split reads were found in $STEP1_UNMASK_BAM" >&2
    fi

    if [[ ! -f "${REFERENCE}.bwt" ]]; then
        echo "Indexing masked / TE-augmented scanning reference..."
        bwa index "$REFERENCE" 2>/dev/null
    fi

    echo "Aligning reads to scanning reference and name-sorting..."
    bwa mem -T 20 -Y -t "$THREADS" "$REFERENCE" "$READ1" "$READ2" | \
        samtools sort -n -@ "$SAMTOOLS_THREADS" -o "$STEP1_BAM" -

    touch "$STEP1_DONE"

    # Alignment changed, so downstream steps must be regenerated.
    rm -f "$STEP2_DONE" "$STEP3_DONE" "$STEP4_DONE"
else
    echo "Step 1/4: Alignment completed, skipping"
fi

# -----------------------------
# Step 2: Fixmate, duplicate removal, name-sort
# -----------------------------
if [[ ! -f "$STEP2_DONE" ]]; then
    echo "Step 2/4: Fixmate and duplicate removal"
    rm -f "$STEP2_DONE"

    if [[ ! -f "$FIXMATE_BAM" ]]; then
        samtools fixmate -@ "$SAMTOOLS_THREADS" -m "$STEP1_BAM" "$FIXMATE_BAM"
    fi

    if [[ ! -f "$POS_SORTED_BAM" ]]; then
        samtools sort -@ "$SAMTOOLS_THREADS" -o "$POS_SORTED_BAM" "$FIXMATE_BAM"
    fi

    if [[ ! -f "$DEDUP_BAM" ]]; then
        samtools rmdup -s "$POS_SORTED_BAM" "$DEDUP_BAM" 2> "${OUT_PREFIX}.rmdup.log"
    fi

    if [[ ! -f "$NAME_SORTED_BAM" ]]; then
        samtools sort -n -@ "$SAMTOOLS_THREADS" -o "$NAME_SORTED_BAM" "$DEDUP_BAM"
    fi

    touch "$STEP2_DONE"

    # Deduplicated/name-sorted BAM changed, so downstream steps must be regenerated.
    rm -f "$STEP3_DONE" "$STEP4_DONE"
else
    echo "Step 2/4: Duplicate removal completed, skipping"
fi

# -----------------------------
# Step 3: Chunk preparation
# -----------------------------
CURRENT_META=$(
    cat <<EOF
CHUNKS=${CHUNKS}
NAME_SORTED_BAM=${NAME_SORTED_BAM}
REFERENCE=${REFERENCE}
READ1=${READ1}
READ2=${READ2}
EOF
)

REBUILD_CHUNKS=0
if [[ ! -f "$STEP3_DONE" ]]; then
    REBUILD_CHUNKS=1
elif [[ ! -f "$STEP3_META" ]]; then
    REBUILD_CHUNKS=1
elif [[ "$CURRENT_META" != "$(cat "$STEP3_META")" ]]; then
    REBUILD_CHUNKS=1
fi

if [[ "$REBUILD_CHUNKS" -eq 1 ]]; then
    echo "Step 3/4: Preparing query-name-sorted BAM chunks"
    rm -f "$STEP3_DONE"
    rm -rf "$CHUNK_DIR"
    mkdir -p "$CHUNK_DIR"

    echo "Extracting mapped read names with MAPQ > 0..."
    samtools view -F 4 -q 1 -@ "$SAMTOOLS_THREADS" "$NAME_SORTED_BAM" | \
        awk '{print $1}' | \
        sort -u > "$READ_NAMES"

    if [[ ! -s "$READ_NAMES" ]]; then
        echo "ERROR: No mapped read names with MAPQ > 0 were found in $NAME_SORTED_BAM" >&2
        exit 1
    fi

    if [[ "$CHUNKS" -gt 1 ]]; then
        echo "Splitting read names into $CHUNKS chunks..."
        split -d -n l/"$CHUNKS" "$READ_NAMES" "${CHUNK_DIR}/chunk_"

        echo "Creating chunk BAMs from the name-sorted deduplicated BAM..."
        for CHUNK_FILE in "${CHUNK_DIR}"/chunk_*; do
            [[ -s "$CHUNK_FILE" ]] || continue

            RAW_CHUNK_BAM="${CHUNK_FILE}.raw.bam"
            FINAL_CHUNK_BAM="${CHUNK_FILE}.bam"

            samtools view -@ "$SAMTOOLS_THREADS" \
                -N "$CHUNK_FILE" \
                -b -o "$RAW_CHUNK_BAM" \
                "$NAME_SORTED_BAM"

            # Explicitly enforce query-name sorting.
            # This is critical because second_pass.py groups adjacent records by read name.
            samtools sort -n -@ "$SAMTOOLS_THREADS" \
                -o "$FINAL_CHUNK_BAM" \
                "$RAW_CHUNK_BAM"

            rm -f "$RAW_CHUNK_BAM"
        done
    else
        ln -sf "$(realpath "$NAME_SORTED_BAM")" "${CHUNK_DIR}/chunk_00.bam"
    fi

    printf "%s\n" "$CURRENT_META" > "$STEP3_META"
    touch "$STEP3_DONE"

    # Chunking changed, so detection must be regenerated.
    rm -f "$STEP4_DONE"
else
    echo "Step 3/4: Chunk preparation completed, skipping"
fi

# -----------------------------
# Step 4: Detection and merging
# -----------------------------
if [[ ! -f "$STEP4_DONE" ]]; then
    echo "Step 4/4: Detecting insertion loci"
    rm -f "$STEP4_DONE"
    rm -rf "${OUT_PREFIX}_secondpass"
    mkdir -p "${OUT_PREFIX}_secondpass"

    echo "Phase 1: Extracting split-read BAM for first_pass.py"

    # Use AWK instead of samtools -N here for broader compatibility and to preserve headers.
    samtools view -h "$DEDUP_BAM" | \
        awk 'NR==FNR {names[$1]; next} /^@/ || ($1 in names)' "$SPLIT_READS" - | \
        samtools view -@ "$SAMTOOLS_THREADS" -b -o "$SPLIT_BAM" -

    # Sort and index the split-read BAM.
    # first_pass.py currently streams BAM records, but this makes the file robust for future fetch-based logic.
    samtools sort -@ "$SAMTOOLS_THREADS" -o "$SPLIT_POS_BAM" "$SPLIT_BAM"
    samtools index "$SPLIT_POS_BAM"

    FIRST_PASS_ARGS=(-b "$SPLIT_POS_BAM" -o "$MERGED_BED")
    if [[ -n "${BED_FILE:-}" ]]; then
        FIRST_PASS_ARGS+=(--bed "$BED_FILE")
    fi

    python3 "${SCRIPT_DIR}/first_pass.py" "${FIRST_PASS_ARGS[@]}"

    echo "Phase 2: Running second_pass.py"

    if [[ "$CHUNKS" -gt 1 ]]; then
        find "$CHUNK_DIR" -name "chunk_*.bam" -print0 | \
            parallel -0 -j "$THREADS" \
                "python3 '${SCRIPT_DIR}/second_pass.py' -b {} -m '$MERGED_BED' -o '${OUT_PREFIX}_secondpass/secondpass_{#}.tsv'"

        python3 "${SCRIPT_DIR}/merge_reports.py" \
            -i "${OUT_PREFIX}_secondpass"/secondpass_*.tsv \
            -o "$STEP4_REPORT"
    else
        python3 "${SCRIPT_DIR}/second_pass.py" \
            -m "$MERGED_BED" \
            -b "$NAME_SORTED_BAM" \
            -o "$STEP4_REPORT"
    fi

    echo "Estimating modal depth from unmasked-reference alignment..."
    samtools depth "$STEP1_UNMASK_BAM" > "${OUT_PREFIX}.depth"

    if [[ ! -s "${OUT_PREFIX}.depth" ]]; then
        echo "ERROR: Depth file is empty: ${OUT_PREFIX}.depth" >&2
        exit 1
    fi

    DEPTH=$(
        shuf -n 1000000 "${OUT_PREFIX}.depth" | \
            cut -f 3 | \
            sort -n | \
            uniq -c | \
            sort -nk1,1 | \
            tail -n1 | \
            awk '{print $2}'
    )

    if [[ -z "$DEPTH" ]]; then
        echo "ERROR: Failed to estimate modal depth." >&2
        exit 1
    fi

    echo "Modal depth estimate: $DEPTH"

    awk -v x="$DEPTH" '
        BEGIN {OFS="\t"}
        NR==1 {
            print
            next
        }
        $5 !~ /Ancestral/ &&
        $5 !~ /Unintact/ &&
        $9 == "Yes" &&
        ($6 + $7 + $8) < 2.5 * x
    ' "$STEP4_REPORT" > "${OUT_PREFIX}.temp.tsv"

    python3 "${SCRIPT_DIR}/id_conversion.py" \
        "${OUT_PREFIX}.temp.tsv" \
        "${SCRIPT_DIR}/IL_CONVERSION.tsv" \
        "$FINAL_TSV"

    touch "$STEP4_DONE"
else
    echo "Step 4/4: Detection completed, skipping"
fi

echo "Pipeline completed."
echo "Final result: $FINAL_TSV"
