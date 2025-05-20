#!/usr/bin/env python3

import sys
import pysam # Re-introducing pysam for efficient VCF reading
import collections

# --- Configuration ---
DEBUG = False # Set to True for detailed debug messages
MAX_SITES_TO_PROCESS = 0 # Set to 0 to process all sites, or a small number for debugging

# --- Helper Functions ---

def print_debug(message):
    """Prints debug messages to stderr if DEBUG is True."""
    if DEBUG:
        print(f"DEBUG: {message}", file=sys.stderr)

def get_genotype_from_sample_data(dp_val, ad_tuple, 
                                  ref_allele_from_final_sites, 
                                  sample_ref_allele, # REF from sample's gVCF (used for AD[0] context)
                                  sample_alt_alleles):
    """
    Applies custom genotyping rules.
    The 'ref_allele_from_final_sites' is used for the indel length comparison AND for reporting homozygous reference.
    'sample_alt_alleles' (from gVCF) are used for reporting homozygous alternate SNVs.

    Returns the genotype as:
    - Reference allele (from final_sites.vcf) if homozygous reference (>=50% ref reads, DP>=4).
    - Sample's Alternate allele (SNV from gVCF) if homozygous alternate SNV (>50% specific SNV ALT reads, DP>=4, SNV_ALT_DP>=4).
    - '.' (missing) if depth < 4, or ambiguous, or if a well-supported indel is found (indel vs ref_from_final_sites).

    Args:
        dp_val (int/None): Total depth (DP) from the sample's VCF record for this site.
        ad_tuple (tuple/None): Allelic depths (AD) from the sample's VCF. AD[0] corresponds to sample_ref_allele.
        ref_allele_from_final_sites (str): REF allele from the final_sites.vcf for indel length comparison and hom-ref reporting.
        sample_ref_allele (str): The reference allele observed in the sample's VCF record (context for AD[0]).
        sample_alt_alleles (tuple/None): Tuple of alternate alleles observed in the sample's VCF record.
    """
    print_debug(f"get_genotype_from_sample_data: Input DP={dp_val}, AD={ad_tuple}, REF_FINAL_SITES='{ref_allele_from_final_sites}', SAMPLE_REF='{sample_ref_allele}', SAMPLE_ALTs={sample_alt_alleles}")

    if dp_val is None or dp_val < 4:
        print_debug(f"get_genotype_from_sample_data: DP {dp_val} < 4 or None. Returning '.'")
        return '.'
    if ad_tuple is None or not ad_tuple: # Checks for None or empty tuple
        print_debug(f"get_genotype_from_sample_data: AD is None or empty. Returning '.'")
        return '.'

    ad_ref = ad_tuple[0] # This AD corresponds to the sample_ref_allele from the gVCF
    ref_freq = ad_ref / dp_val
    print_debug(f"get_genotype_from_sample_data: Calculated REF_FREQ (for sample_ref_allele '{sample_ref_allele}') = {ref_freq:.4f} (AD_REF={ad_ref}/DP={dp_val})")

    # Rule: Homozygous reference (output ref_allele_from_final_sites) if >=50% reads support the sample_ref_allele
    if ref_freq >= 0.90:
        print_debug(f"get_genotype_from_sample_data: REF_FREQ >= 0.90. Returning REF_FROM_FINAL_SITES: '{ref_allele_from_final_sites}'")
        return ref_allele_from_final_sites # MODIFIED: Report REF from final_sites.vcf

    # If not homozygous reference, check alternate alleles from the sample's gVCF.
    if sample_alt_alleles:
        print_debug(f"get_genotype_from_sample_data: Checking SAMPLE_ALTs: {sample_alt_alleles}")
        for idx, current_sample_alt_allele in enumerate(sample_alt_alleles):
            if current_sample_alt_allele is None or current_sample_alt_allele == '.':
                print_debug(f"get_genotype_from_sample_data: Skipping sample ALT '.' or None at index {idx}")
                continue

            ad_current_alt = 0
            ad_index = idx + 1 # Index in the AD tuple (0 is ref)
            if len(ad_tuple) > ad_index:
                ad_val_from_tuple = ad_tuple[ad_index]
                if isinstance(ad_val_from_tuple, int):
                    ad_current_alt = ad_val_from_tuple
                else:
                     print_debug(f"get_genotype_from_sample_data: AD value for sample ALT '{current_sample_alt_allele}' at index {ad_index} is not an int: {ad_val_from_tuple}. Treating as 0.")
            else:
                 print_debug(f"get_genotype_from_sample_data: AD tuple {ad_tuple} too short for sample ALT '{current_sample_alt_allele}' at index {ad_index}. Treating as 0.")
            
            print_debug(f"get_genotype_from_sample_data: Checking sample ALT '{current_sample_alt_allele}': AD_ALT={ad_current_alt}")

            # Check if this alternate allele from the sample is well-supported
            if ref_freq<0.90 and ref_freq>0.10:
                return '.'
            if ad_current_alt >= 4:
                current_alt_freq = ad_current_alt / dp_val
                print_debug(f"get_genotype_from_sample_data: Sample ALT '{current_sample_alt_allele}': AD_ALT >= 4. Calculated ALT_FREQ = {current_alt_freq:.4f}")
                if current_alt_freq >= 0.90:
                    # This sample ALT is homozygous by frequency and depth.
                    # Now check if it's an indel with respect to the REF from final_sites.vcf
                    is_indel = (len(ref_allele_from_final_sites) != len(current_sample_alt_allele))
                    print_debug(f"get_genotype_from_sample_data: Sample ALT '{current_sample_alt_allele}': REF_FINAL_SITES_LEN={len(ref_allele_from_final_sites)}, CURRENT_SAMPLE_ALT_LEN={len(current_sample_alt_allele)}, IsIndel={is_indel}")
                    
                    if is_indel:
                        print_debug(f"get_genotype_from_sample_data: Well-supported indel '{current_sample_alt_allele}' (vs REF_FINAL_SITES) found. Returning '.'")
                        return '.' # Report as missing if it's a well-supported indel
                    else:
                        # It's a well-supported SNV
                        print_debug(f"get_genotype_from_sample_data: Well-supported SNV. Returning '{current_sample_alt_allele}'")
                        return current_sample_alt_allele # Report the SNV alternate allele from sample
                else:
                    print_debug(f"get_genotype_from_sample_data: Sample ALT '{current_sample_alt_allele}': ALT_FREQ < 0.90.")
            else:
                 print_debug(f"get_genotype_from_sample_data: Sample ALT '{current_sample_alt_allele}': AD_ALT < 4.")

    print_debug(f"get_genotype_from_sample_data: No allele met criteria or ambiguous. Returning '.'")
    return '.'


def main(final_sites_vcf_path, sample_gvl_vcf_path, sample_id):
    """
    Uses pysam to read target coordinates and fetch data from the sample's gVCF.
    Outputs one genotype per line for the given sample.
    """
    print_debug(f"--- Starting Genotyping for Sample: {sample_id} (pysam version) ---")
    print_debug(f"Sites VCF: {final_sites_vcf_path}")
    print_debug(f"Sample gVCF: {sample_gvl_vcf_path}")

    # Stores dicts: {'chrom': str, 'pos': int, 'ref_from_final_sites': str, 'alt_from_final_sites': str, (optional)'genotype': '.'}
    ordered_target_sites_info = []
    try:
        print_debug(f"Reading target coordinates and REF from {final_sites_vcf_path} using pysam...")
        sites_vcf_pysam = pysam.VariantFile(final_sites_vcf_path, "r")
        sites_processed_count = 0
        for site_record in sites_vcf_pysam.fetch():
            if MAX_SITES_TO_PROCESS > 0 and sites_processed_count >= MAX_SITES_TO_PROCESS:
                print_debug(f"Reached MAX_SITES_TO_PROCESS ({MAX_SITES_TO_PROCESS}). Stopping coordinate reading.")
                break
            
            # final_sites.vcf should contain biallelic SNVs
            if not site_record.alts or len(site_record.alts) != 1 or site_record.alts[0] is None or site_record.alts[0] == '.':
                print_debug(f"Warning: Skipping site {site_record.chrom}:{site_record.pos} in {final_sites_vcf_path} - not a valid biallelic SNV (REF: {site_record.ref}, ALT: {site_record.alts}).")
                ordered_target_sites_info.append({'chrom': site_record.chrom, 'pos': site_record.pos, 'genotype': '.'}) # Mark for printing '.'
                continue
            
            ordered_target_sites_info.append({
                'chrom': site_record.chrom, 
                'pos': site_record.pos,
                'ref_from_final_sites': site_record.ref, # REF from final_sites.vcf
                'alt_from_final_sites': site_record.alts[0]  # ALT from final_sites.vcf (for context)
            })
            sites_processed_count += 1
        print_debug(f"Read {len(ordered_target_sites_info)} target coordinates.")
    except Exception as e:
        print_debug(f"FATAL: Error reading coordinates from final_sites.vcf '{final_sites_vcf_path}' with pysam: {e}")
        return
    finally:
        if 'sites_vcf_pysam' in locals() and sites_vcf_pysam:
            sites_vcf_pysam.close()

    if not ordered_target_sites_info:
        print_debug(f"Warning: No target coordinates to process from {final_sites_vcf_path}")
        return

    try:
        sample_vcf_pysam = pysam.VariantFile(sample_gvl_vcf_path, "r")
    except Exception as e:
        print_debug(f"FATAL: Could not open sample gVL VCF '{sample_gvl_vcf_path}' with pysam for sample '{sample_id}'. Error: {e}")
        for _ in ordered_target_sites_info: print('.')
        return

    print_debug(f"--- Determining Genotypes for {len(ordered_target_sites_info)} Target Coordinates ---")
    for target_site_spec in ordered_target_sites_info:
        genotype_to_print = '.'
        
        if 'genotype' in target_site_spec and target_site_spec['genotype'] == '.':
            print(genotype_to_print)
            print_debug(f"Target coordinate {target_site_spec['chrom']}:{target_site_spec['pos']} was pre-marked as invalid. Genotype: '.'")
            print_debug("-" * 20)
            continue
            
        chrom = target_site_spec['chrom']
        pos = target_site_spec['pos']
        ref_allele_from_final_sites = target_site_spec['ref_from_final_sites'] # This is the key change for this version

        print_debug(f"Processing target coordinate: ({chrom}, {pos}), REF_from_final_sites: '{ref_allele_from_final_sites}'")

        try:
            fetched_records = list(sample_vcf_pysam.fetch(chrom, pos - 1, pos))
            if fetched_records:
                sample_site_record = fetched_records[0]
                
                if sample_id not in sample_site_record.samples:
                    print_debug(f"Warning: Sample ID '{sample_id}' not found in pysam record for {chrom}:{pos}. Defaulting to '.'")
                    print(genotype_to_print)
                    print_debug("-" * 20)
                    continue

                fmt = sample_site_record.samples[sample_id]
                dp_val = fmt.get('DP')
                ad_from_pysam = fmt.get('AD')

                processed_ad_tuple = None
                if isinstance(ad_from_pysam, int):
                    processed_ad_tuple = (ad_from_pysam,)
                elif isinstance(ad_from_pysam, tuple):
                    processed_ad_tuple = tuple(ad if isinstance(ad, int) else 0 for ad in ad_from_pysam)
                
                sample_ref_from_gvl = sample_site_record.ref
                sample_alts_from_gvl = sample_site_record.alts

                print_debug(f"Data found for ({chrom}, {pos}) in gVCF: DP={dp_val}, AD(raw)={ad_from_pysam} -> AD(proc)={processed_ad_tuple}, SAMPLE_REF='{sample_ref_from_gvl}', SAMPLE_ALTs={sample_alts_from_gvl}")

                # Pass REF from final_sites.vcf for indel check AND hom-ref reporting
                genotype_to_print = get_genotype_from_sample_data(
                    dp_val=dp_val,
                    ad_tuple=processed_ad_tuple,
                    ref_allele_from_final_sites=ref_allele_from_final_sites, # Passed for hom-ref reporting
                    sample_ref_allele=sample_ref_from_gvl, # Passed for AD[0] context
                    sample_alt_alleles=sample_alts_from_gvl
                )
            else:
                print_debug(f"Coordinate ({chrom}, {pos}) not found in sample gVCF using pysam.fetch. Genotype: '.'")
        except Exception as e:
            print_debug(f"Error during pysam fetch or processing for {chrom}:{pos} in {sample_gvl_vcf_path}: {e}. Defaulting to '.'")

        print(genotype_to_print)
        print_debug(f"Final genotype for ({chrom}, {pos}): {genotype_to_print}")
        print_debug("-" * 20)
    
    if 'sample_vcf_pysam' in locals() and sample_vcf_pysam:
        sample_vcf_pysam.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python genotype_single_sample.py <final_sites.vcf.gz> <sample.gvl.vcf.gz> <sample_id>", file=sys.stderr)
        sys.exit(1)
    final_sites_vcf_path = sys.argv[1]
    sample_gvl_vcf_path = sys.argv[2]
    sample_id_arg = sys.argv[3]
    if MAX_SITES_TO_PROCESS > 0:
        print_debug(f"WARNING: Script will only process the first {MAX_SITES_TO_PROCESS} sites for debugging.")
    main(final_sites_vcf_path, sample_gvl_vcf_path, sample_id_arg)

