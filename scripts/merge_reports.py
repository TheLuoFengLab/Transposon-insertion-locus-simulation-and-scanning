#!/usr/bin/env python3

import pandas as pd
import argparse

def process_reports(input_files):
    """Aggregate counts without coordinate merging"""
    # Read all files
    dfs = []
    for f in input_files:
        df = pd.read_csv(f, sep='\t', dtype={
            'Chrom': str, 'Start': int, 'End': int,
            'TE_Family': str, 'Type': str,
            'TypeA': int, 'TypeB': int, 'TypeC': int, 'PASS_FILTER': str
        })
        dfs.append(df)
    
    # Combine all data
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Group by exact coordinates and TE family
    grouped = combined_df.groupby(['Chrom', 'Start', 'End', 'TE_Family', 'Type'])

    
    merged_df = grouped.agg({
        'TypeA': 'sum',
        'TypeB': 'sum',
        'TypeC': 'sum'
    }).reset_index()
    
    # Calculate PASS_FILTER after aggregation
    merged_df['PASS_FILTER'] = (merged_df['TypeB'] + merged_df['TypeA'] > 1).map({True: 'Yes', False: 'No'})
    
    return merged_df[['Chrom', 'Start', 'End', 'TE_Family', 'Type',
                     'TypeA', 'TypeB', 'TypeC', 'PASS_FILTER']]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs='+', required=True, help='Input TSV files')
    parser.add_argument('-o', required=True, help='Output file')
    args = parser.parse_args()

    result = process_reports(args.i)
    result.to_csv(args.o, sep='\t', header=True, index=False)