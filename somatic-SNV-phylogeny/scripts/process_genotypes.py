#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

def print_stderr(message):
    """Prints messages to stderr."""
    print(message, file=sys.stderr)

def filter_genotypes(input_tsv_path, output_filtered_tsv_path,
                     max_variant_missing_pct=0.10,
                     max_sample_missing_pct=0.10):
    """
    Filters variants and samples based on missing data percentages,
    and then removes invariant sites.

    Args:
        input_tsv_path (str): Path to the input concatenated genotype TSV file.
        output_filtered_tsv_path (str): Path to save the filtered TSV file.
        max_variant_missing_pct (float): Maximum allowed missing data percentage for a variant.
        max_sample_missing_pct (float): Maximum allowed missing data percentage for a sample.

    Returns:
        pandas.DataFrame: The filtered DataFrame, or None if an error occurs.
    """
    print_stderr(f"Loading genotype data from: {input_tsv_path}")
    try:
        # Read the TSV file. Assuming the first line is the header.
        df = pd.read_csv(input_tsv_path, sep='\t', header=0)
        # The user indicated the #CHROM renaming was commented out in their version.
        # if df.columns[0].startswith('#'):
        #     df.rename(columns={df.columns[0]: df.columns[0][1:]}, inplace=True)
        
        print_stderr(f"Initial data shape: {df.shape[0]} variants, {df.shape[1]-4} samples (excluding CHROM, POS, REF, ALT)")

    except FileNotFoundError:
        print_stderr(f"Error: Input TSV file not found at {input_tsv_path}")
        return None
    except Exception as e:
        print_stderr(f"Error loading TSV file: {e}")
        return None

    if df.shape[1] < 5: # CHROM, POS, REF, ALT, Sample1
        print_stderr("Error: Input TSV must have at least 5 columns (CHROM, POS, REF, ALT, Sample1...).")
        return None

    site_info_cols = ['CHROM', 'POS', 'REF', 'ALT']
    for col in site_info_cols:
        if col not in df.columns:
            print_stderr(f"Error: Expected column '{col}' not found in the input TSV.")
            print_stderr(f"Found columns: {df.columns.tolist()}")
            return None
            
    genotype_cols = df.columns[4:].tolist()
    num_initial_samples = len(genotype_cols)
    if num_initial_samples == 0:
        print_stderr("Error: No sample genotype columns found.")
        return None

    print_stderr(f"Identified {num_initial_samples} initial sample columns: {genotype_cols}")

    # --- 1. Filter variants (rows) with > max_variant_missing_pct missing data ---
    print_stderr(f"\nFiltering variants with > {max_variant_missing_pct*100:.1f}% missing data...")
    # Ensure genotype columns are treated as strings to correctly identify '.'
    for col in genotype_cols:
        df[col] = df[col].astype(str)
        
    df['missing_count_variant'] = df[genotype_cols].apply(lambda row: (row == '.').sum(), axis=1)
    df['missing_pct_variant'] = df['missing_count_variant'] / num_initial_samples
    
    variants_before_filter1 = df.shape[0]
    df_after_variant_missing_filter = df[df['missing_pct_variant'] <= max_variant_missing_pct].copy()
    variants_after_filter1 = df_after_variant_missing_filter.shape[0]
    print_stderr(f"Variants before missingness filter: {variants_before_filter1}")
    print_stderr(f"Variants after missingness filter: {variants_after_filter1} (Removed {variants_before_filter1 - variants_after_filter1})")

    if df_after_variant_missing_filter.empty:
        print_stderr("Error: No variants remained after variant missingness filtering. Cannot proceed.")
        return None
    df_after_variant_missing_filter.drop(columns=['missing_count_variant', 'missing_pct_variant'], inplace=True)

    # --- 2. Filter samples (columns) with > max_sample_missing_pct missing data ---
    print_stderr(f"\nFiltering samples with > {max_sample_missing_pct*100:.1f}% missing data (on {variants_after_filter1} variants)...")
    
    missing_per_sample = df_after_variant_missing_filter[genotype_cols].apply(lambda col: (col == '.').sum())
    missing_pct_per_sample = missing_per_sample / variants_after_filter1 

    samples_to_keep = missing_pct_per_sample[missing_pct_per_sample <= max_sample_missing_pct].index.tolist()
    samples_removed = [s for s in genotype_cols if s not in samples_to_keep]

    print_stderr(f"Samples before missingness filter: {num_initial_samples}")
    print_stderr(f"Samples after missingness filter: {len(samples_to_keep)} (Removed {len(samples_removed)}: {samples_removed})")

    if not samples_to_keep:
        print_stderr("Error: No samples remained after sample missingness filtering. Cannot proceed.")
        return None

    current_genotype_cols = samples_to_keep
    df_after_sample_missing_filter = df_after_variant_missing_filter[site_info_cols + current_genotype_cols].copy()
    
    print_stderr(f"Data shape after sample missingness filter: {df_after_sample_missing_filter.shape[0]} variants, {len(current_genotype_cols)} samples.")

    if df_after_sample_missing_filter.empty:
        print_stderr("Error: DataFrame is empty after sample missingness filtering. Cannot proceed.")
        return None
        
    # --- 3. Filter invariant sites ---
    print_stderr(f"\nFiltering invariant sites (based on {len(current_genotype_cols)} samples)...")
    
    def is_invariant(row):
        genotypes_for_current_samples = row[current_genotype_cols]
        # Ensure comparison is with string '.'
        unique_non_missing_alleles = set(g for g in genotypes_for_current_samples.astype(str) if g != '.')
        return len(unique_non_missing_alleles) < 2 

    invariant_mask = df_after_sample_missing_filter.apply(is_invariant, axis=1)
    df_final_filtered = df_after_sample_missing_filter[~invariant_mask].copy()

    variants_before_invariant_filter = df_after_sample_missing_filter.shape[0]
    variants_after_invariant_filter = df_final_filtered.shape[0]
    print_stderr(f"Variants before invariant site filter: {variants_before_invariant_filter}")
    print_stderr(f"Variants after invariant site filter: {variants_after_invariant_filter} (Removed {variants_before_invariant_filter - variants_after_invariant_filter})")

    if df_final_filtered.empty:
        print_stderr("Error: No variants remained after invariant site filtering. Cannot proceed.")
        return None

    print_stderr(f"\nFinal filtered data shape: {df_final_filtered.shape[0]} variants, {len(current_genotype_cols)} samples.")
    
    try:
        df_final_filtered.to_csv(output_filtered_tsv_path, sep='\t', index=False)
        print_stderr(f"Filtered TSV saved to: {output_filtered_tsv_path}")
    except Exception as e:
        print_stderr(f"Error saving filtered TSV: {e}")
        return None
        
    return df_final_filtered


def convert_to_phylip(filtered_df, output_phylip_path):
    """
    Converts the filtered genotype DataFrame to PHYLIP sequential format.
    Replaces '.' with '?' for missing data.

    Args:
        filtered_df (pandas.DataFrame): The filtered DataFrame.
        output_phylip_path (str): Path to save the PHYLIP file.
    """
    if filtered_df is None or filtered_df.empty:
        print_stderr("Cannot convert to PHYLIP: input DataFrame is empty or None.")
        return

    site_info_cols = ['CHROM', 'POS', 'REF', 'ALT']
    sample_names = filtered_df.columns[len(site_info_cols):].tolist()
    num_taxa = len(sample_names)
    num_sites = filtered_df.shape[0]

    if num_taxa == 0 or num_sites == 0:
        print_stderr("Cannot convert to PHYLIP: No samples or no sites after filtering.")
        return

    print_stderr(f"\nConverting to PHYLIP format: {num_taxa} taxa, {num_sites} sites. Missing data '.' will be converted to '?'.")

    # Create a copy for modification to avoid changing the DataFrame passed to it
    df_for_phylip = filtered_df.copy()

    # Replace '.' with '?' in genotype columns for PHYLIP output
    for col in sample_names:
        df_for_phylip[col] = df_for_phylip[col].astype(str).replace({'.': '?'})

    with open(output_phylip_path, 'w') as f_phy:
        f_phy.write(f" {num_taxa} {num_sites}\n")

        for sample_name in sample_names:
            if len(sample_name) > 10:
                padded_name = sample_name[:10]
                print_stderr(f"Warning: Sample name '{sample_name}' truncated to '{padded_name}' for PHYLIP format.")
            else:
                padded_name = sample_name.ljust(10)

            # Sequence is now composed of A, C, G, T, or ?
            sequence = "".join(df_for_phylip[sample_name].tolist()) # Already strings
            
            f_phy.write(f"{padded_name}{sequence}\n")

    print_stderr(f"PHYLIP file saved to: {output_phylip_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter a concatenated genotype TSV file, remove invariant sites, and convert it to PHYLIP format."
    )
    parser.add_argument(
        "input_tsv",
        help="Path to the input concatenated genotype TSV file (e.g., final_genotypes.tsv)."
    )
    parser.add_argument(
        "-o", "--output_prefix",
        default="filtered_genotypes",
        help="Prefix for the output files (e.g., 'filtered_genotypes' will create "
             "'filtered_genotypes.tsv' and 'filtered_genotypes.phy'). Default: 'filtered_genotypes'."
    )
    parser.add_argument(
        "--max_variant_missing", type=float, default=0.10,
        help="Maximum allowed missing data percentage for a variant (0.0 to 1.0). Default: 0.10 (10%%)."
    )
    parser.add_argument(
        "--max_sample_missing", type=float, default=0.10,
        help="Maximum allowed missing data percentage for a sample (0.0 to 1.0). Default: 0.10 (10%%)."
    )

    args = parser.parse_args()

    output_filtered_tsv = f"{args.output_prefix}.filtered.tsv"
    output_phylip = f"{args.output_prefix}.phy"

    filtered_df = filter_genotypes(
        args.input_tsv,
        output_filtered_tsv, 
        args.max_variant_missing,
        args.max_sample_missing
    )

    if filtered_df is not None and not filtered_df.empty:
        convert_to_phylip(filtered_df, output_phylip)
    else:
        print_stderr("Skipping PHYLIP conversion due to errors or no data after filtering.")

if __name__ == "__main__":
    main()

