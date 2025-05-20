#!/bin/bash

# Exit script on error and treat unset variables as errors
set -e -u -o pipefail

# --- Configuration ---
PYTHON_CMD="python3" # Command to run python
FILTER_SCRIPT_PATH="scripts/process_genotypes.py" # Path to your Python filtering script
BCFTOOLS_CMD="bcftools" # Path to bcftools, if not in PATH

# --- Usage ---
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <final_sites_vcf> <input_genotype_cols_directory> <output_prefix>"
    echo "  <final_sites_vcf>: Path to the final_sites.vcf.gz from the WGS pipeline (contains CHROM,POS,REF,ALT)."
    echo "  <input_genotype_cols_directory>: Directory containing individual sample genotype column files (e.g., output of genotype_single_sample.py)."
    echo "  <output_prefix>: Prefix for output files (e.g., 'my_analysis')."
    echo "  <missing_variant>: ratio of missing genotypes allowed for a variant (float, e.g., '0.05')"
    echo "  <missing_sample>: ratio of missing variant genotypes allowed for a sample (float, e.g., '0.05')"
    exit 1
fi

FINAL_SITES_VCF="$1"
INPUT_GENO_COLS_DIR="$2" # Renamed for clarity
OUTPUT_PREFIX="$3"
MISSING_VARIANT=$4
MISSING_SAMPLE=$5

# --- Output Directories and Files ---
CONCAT_DIR="01_concatenated"
FILTER_DIR="02_filtered_converted"
LOG_DIR="logs"

CONCAT_TSV="${CONCAT_DIR}/${OUTPUT_PREFIX}.concatenated.tsv"
FILTERED_TSV_PREFIX="${FILTER_DIR}/${OUTPUT_PREFIX}" # Python script appends .filtered.tsv and .phy
FINAL_PHY_FILE="${FILTERED_TSV_PREFIX}.phy" # Path to the final PHYLIP file

mkdir -p "$CONCAT_DIR" "$FILTER_DIR" "$LOG_DIR"

# --- Sanity Checks ---
if [ ! -f "$FINAL_SITES_VCF" ]; then
    echo "Error: Final sites VCF file '$FINAL_SITES_VCF' not found."
    exit 1
fi
if [ ! -d "$INPUT_GENO_COLS_DIR" ]; then
    echo "Error: Input genotype columns directory '$INPUT_GENO_COLS_DIR' not found."
    exit 1
fi
if [ ! -f "$FILTER_SCRIPT_PATH" ]; then
    echo "Error: Python filter script '$FILTER_SCRIPT_PATH' not found."
    exit 1
fi
if ! command -v $PYTHON_CMD &> /dev/null; then
    echo "Error: Python command '$PYTHON_CMD' not found."
    exit 1
fi
if ! command -v $BCFTOOLS_CMD &> /dev/null; then
    echo "Error: bcftools command '$BCFTOOLS_CMD' not found."
    exit 1
fi
if ! command -v awk &> /dev/null; then echo "Error: awk not found."; exit 1; fi
if ! command -v paste &> /dev/null; then echo "Error: paste not found."; exit 1; fi
if ! command -v basename &> /dev/null; then echo "Error: basename not found."; exit 1; fi
if ! command -v head &> /dev/null; then echo "Error: head not found."; exit 1; fi
if ! command -v wc &> /dev/null; then echo "Error: wc not found."; exit 1; fi


echo "INFO: Checks passed. Final sites VCF: $FINAL_SITES_VCF, Input genotype columns directory: $INPUT_GENO_COLS_DIR, Output prefix: $OUTPUT_PREFIX"

# --- Step 1: Concatenate site information and individual sample genotype column files ---
echo "INFO: Step 1 - Concatenating site info and sample genotype columns..."

# Create a temporary directory for intermediate files
TMP_CONCAT_DIR=$(mktemp -d -p "$CONCAT_DIR" concat_tmp.XXXXXX)
trap 'rm -rf "$TMP_CONCAT_DIR"' EXIT # Ensure cleanup on script exit

# Extract site information (CHROM, POS, REF, ALT) from the FINAL_SITES_VCF
# This VCF should already be filtered for biallelic SNVs.
echo "INFO: Extracting site information from '$FINAL_SITES_VCF'..."
$BCFTOOLS_CMD query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$FINAL_SITES_VCF" > "${TMP_CONCAT_DIR}/site_info.tmp"

if [ ! -s "${TMP_CONCAT_DIR}/site_info.tmp" ]; then
    echo "Error: Failed to extract site information from '$FINAL_SITES_VCF', or it contained no variant data."
    exit 1
fi
num_sites_from_vcf=$(wc -l < "${TMP_CONCAT_DIR}/site_info.tmp")
echo "INFO: Extracted $num_sites_from_vcf sites."


# Prepare header and collect genotype column files
genotype_files_to_paste=() # Array to store paths to genotype column files
sample_names_for_header=() # Array to store sample names for the header

# Find sample genotype column files (assuming they are .tsv, adjust if different)
# The output of genotype_single_sample.py was named ${SAMPLE_ID}.geno.txt in the previous pipeline.
# Let's assume the files in INPUT_GENO_COLS_DIR are named like <sample_id>.tsv or <sample_id>.geno.txt
shopt -s nullglob # Important for empty directory case
sample_col_files=("$INPUT_GENO_COLS_DIR"/*.tsv "$INPUT_GENO_COLS_DIR"/*.txt)
shopt -u nullglob

if [ ${#sample_col_files[@]} -eq 0 ]; then
    echo "Error: No .tsv or .txt files found in '$INPUT_GENO_COLS_DIR'."
    exit 1
fi
echo "INFO: Found ${#sample_col_files[@]} potential sample genotype column files."

for sample_col_file in "${sample_col_files[@]}"; do
    if [ ! -f "$sample_col_file" ]; then
        echo "Warning: File '$sample_col_file' not found (this shouldn't happen). Skipping."
        continue
    fi
    # Infer sample name from filename (remove .tsv or .geno.txt suffix)
    sample_name_with_ext=$(basename "$sample_col_file")
    sample_name=${sample_name_with_ext%.tsv}
    sample_name=${sample_name%.geno.txt} # Handle both potential extensions

    echo "INFO: Processing sample: $sample_name from file $sample_col_file"
    
    # Check line count consistency
    num_lines_in_sample_file=$(wc -l < "$sample_col_file")
    if [ "$num_lines_in_sample_file" -ne "$num_sites_from_vcf" ]; then
        echo "Error: Line count mismatch for sample '$sample_name'."
        echo "       Expected $num_sites_from_vcf lines (from $FINAL_SITES_VCF), but '$sample_col_file' has $num_lines_in_sample_file lines."
        echo "       Ensure each sample genotype file corresponds exactly to the sites in $FINAL_SITES_VCF."
        exit 1
    fi
    
    sample_names_for_header+=("$sample_name")
    genotype_files_to_paste+=("$sample_col_file") # Add the original file path directly
done

# Construct the header line
header_line="CHROM\tPOS\tREF\tALT"
for s_name in "${sample_names_for_header[@]}"; do
    header_line="${header_line}\t${s_name}"
done

# Write the constructed header to the concatenated TSV file
echo -e "$header_line" > "$CONCAT_TSV"

# Paste the site information and all genotype columns together
# The order of genotype_files_to_paste matches the order of sample_col_files array
echo "INFO: Pasting site information and genotype columns..."
paste "${TMP_CONCAT_DIR}/site_info.tmp" "${genotype_files_to_paste[@]}" >> "$CONCAT_TSV"

scripts/id_conversion.sh $CONCAT_TSV

# Check if concatenated file was created and has content beyond header
if [ ! -s "$CONCAT_TSV" ] || [ $(wc -l < "$CONCAT_TSV") -le 1 ]; then
    echo "Error: Concatenated TSV file '$CONCAT_TSV' is empty or contains only the header."
    exit 1
fi

echo "INFO: Concatenation complete. Output: $CONCAT_TSV"
echo "INFO: Number of lines in concatenated TSV (incl. header): $(wc -l < "$CONCAT_TSV")"
echo "INFO: Number of columns in concatenated TSV: $(head -n1 "$CONCAT_TSV" | awk -F'\t' '{print NF}')"


# --- Step 2: Filter concatenated TSV and convert to PHYLIP ---
echo -e "\nINFO: Step 2 - Filtering genotypes and converting to PHYLIP format..."
echo "INFO: Using Python script: $FILTER_SCRIPT_PATH"
echo "INFO: Input for Python script: $CONCAT_TSV"
echo "INFO: Output prefix for Python script: $FILTERED_TSV_PREFIX"

# Assuming default missing data thresholds from the python script (10%)
# You can add arguments here if you want to make them configurable in the bash script
$PYTHON_CMD "$FILTER_SCRIPT_PATH" "$CONCAT_TSV" -o "$FILTERED_TSV_PREFIX" --max_variant_missing $MISSING_VARIANT --max_sample_missing $MISSING_SAMPLE

if [ ! -f "$FINAL_PHY_FILE" ]; then
    echo "Error: PHYLIP file '$FINAL_PHY_FILE' was not created by the Python script."
    exit 1
fi
if [ ! -s "$FINAL_PHY_FILE" ]; then
    echo "Error: PHYLIP file '$FINAL_PHY_FILE' is empty."
    exit 1
fi

echo "INFO: Filtering and PHYLIP conversion complete. Output PHYLIP: $FINAL_PHY_FILE"

# --- Step 3: Suggest IQ-TREE2 command ---
echo -e "\nINFO: Step 3 - IQ-TREE2 command suggestion"
echo "---------------------------------------------------------------------"
echo "Suggested IQ-TREE2 command for SNP data (using ModelFinder Plus):"
echo ""
echo "iqtree2 -s \"$FINAL_PHY_FILE\" \\"
echo "        -m MFP \\"
echo "        -T AUTO \\"
echo "        --prefix \"${OUTPUT_PREFIX}_iqtree_mfp\" \\"
echo "        -B 1000"
echo ""
echo "Alternatively, for a common fixed model for SNPs (GTR+ASC+G):"
iqtree2 -s \"$FINAL_PHY_FILE\" \
        -m GTR+ASC+G \
        -T AUTO \
        --prefix \"${OUTPUT_PREFIX}_iqtree_gtr_asc_g\" \
        -B 1000
echo ""
echo "Explanation of parameters:"
echo "  -s: Input alignment file in PHYLIP format."
echo "  -m MFP: ModelFinder Plus to automatically select the best-fit model and apply ascertainment bias correction (ASC) if appropriate for SNP data."
echo "  -m GTR+ASC+G: General Time Reversible model, with ascertainment bias correction for SNPs, and gamma-distributed rate heterogeneity."
echo "  -T AUTO: Automatically determine the optimal number of threads."
echo "  --prefix: Prefix for output files from IQ-TREE2."
echo "  -B 1000: Perform 1000 standard bootstrap replicates."
echo "  (Consider using -alrt 1000 for SH-aLRT branch support if bootstrapping is too slow)."
echo "---------------------------------------------------------------------"

echo -e "\nINFO: Pipeline finished successfully!"

# Cleanup temporary directory (trap will also do this on exit)
rm -rf "$TMP_CONCAT_DIR"

