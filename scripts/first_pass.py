#!/usr/bin/env python3

import pysam
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='First pass IL detection')
    parser.add_argument('-b', '--bam', required=True, help='Input BAM file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    parser.add_argument('--bed', help='Optional BED file with known ILs')
    return parser.parse_args()

def is_te(contig):
    return contig.startswith('TE')

def get_te_family(te_id):
    return te_id.split('.')[0] if '.' in te_id else te_id

class FirstPassProcessor:
    def __init__(self, bam_path, bed_file=None):
        self.bam = pysam.AlignmentFile(bam_path, 'rb')
        self.te_contigs = {ctg for ctg in self.bam.references if is_te(ctg)}
        self.merged_ils = defaultdict(list)
        self.bed_regions = self.load_bed(bed_file) if bed_file else None

    def load_bed(self, bed_file):
        regions = defaultdict(list)
        with open(bed_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 5: continue
                chrom, start, end, family_id, bed_type = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
                regions[chrom].append((max(0, start-20), end+20, family_id, bed_type))
        return regions

    def calculate_split_position(self, read, sa_entry):
        """Calculate precise split position with clipping handling"""
        # Parse SA entry
        sa_fields = sa_entry.split(',')
        if len(sa_fields) < 6:
            return None
        
        try:
            sa_contig = sa_fields[0]
            sa_pos = int(sa_fields[1])  # 1-based POS
            sa_strand = sa_fields[2]
            sa_cigar_str = sa_fields[3]
            sa_mapq = int(sa_fields[4])
        except (ValueError, IndexError):
            return None

        # Check MAPQ requirements
        if read.mapping_quality <= 0 or sa_mapq <= 0:
            return None

        # Determine TE status
        primary_is_te = is_te(read.reference_name)
        sa_is_te = is_te(sa_contig)
        if primary_is_te == sa_is_te:
            return None  # Not a TE/non-TE split

        # Parse CIGAR strings
        def parse_cigar(cigar_str):
            """Parse CIGAR string into list of (operation, length) tuples"""
            import re
            ops = []
            op_map = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4, 
                    'H':5, 'P':6, '=':7, 'X':8}
            for match in re.finditer(r'(\d+)([MIDNSHPX=])', cigar_str):
                length, op = match.groups()
                ops.append((op_map[op], int(length)))
            return ops

        try:
            sa_cigar = parse_cigar(sa_cigar_str)
            primary_cigar = read.cigartuples or []
        except (KeyError, ValueError):
            return None

        # Identify non-TE component and its CIGAR
        if primary_is_te:
            non_te_cigar = sa_cigar
            # Calculate SA alignment coordinates
            sa_ref_length = sum(length for op, length in sa_cigar if op in (0,2,3,7,8))
            non_te_start = sa_pos - 1  # Convert to 0-based
            non_te_end = non_te_start + sa_ref_length
        else:
            non_te_cigar = primary_cigar
            non_te_start = read.reference_start
            non_te_end = read.reference_end

        # Analyze clipping in non-TE component
        def get_clipping(cigar):
            """Return (leading_clip, trailing_clip) in bases"""
            leading = trailing = 0
            if cigar:
                first_op, first_len = cigar[0]
                if first_op in (4,5):  # S/H
                    leading = first_len
                last_op, last_len = cigar[-1]
                if last_op in (4,5):  # S/H
                    trailing = last_len
            return (leading, trailing)

        lead_clip, trail_clip = get_clipping(non_te_cigar)
        
        # Handle excessive clipping
        if lead_clip >= 15 and trail_clip >= 15:
            return None
        
        # Determine split position based on clipping
        if lead_clip >= 15:
            return non_te_start  # Start of non-TE alignment
        elif trail_clip >= 15:
            return non_te_end - 1  # End of non-TE alignment
        
        # Original strand-based calculation
        if primary_is_te:
            # TE -> non-TE split
            if read.is_reverse:
                return read.reference_start  # Start of TE alignment
            else:
                return read.reference_end - 1  # End of TE alignment
        else:
            # non-TE -> TE split
            if sa_strand == '-':
                # TE alignment end (1-based)
                sa_ref_length = sum(length for op, length in sa_cigar if op in (0,2,3,7,8))
                return sa_pos - 1 + sa_ref_length - 1  # 0-based end
            else:
                return sa_pos - 1  # 0-based start of TE alignment

    def process_reads(self):
        il_candidates = defaultdict(lambda: defaultdict(list))
        
        # Original first_pass read processing
        for read in self.bam:
            if read.is_secondary or read.is_supplementary or read.mapping_quality == 0:
                continue
            
            if read.has_tag('SA'):
                for sa_entry in read.get_tag('SA').split(';'):
                    if not sa_entry: continue
                    sa_fields = sa_entry.split(',')
                    if len(sa_fields) < 6: continue
                    sa_contig, sa_pos = sa_fields[0], int(sa_fields[1])
                    sa_mapq = int(sa_fields[4])
                    if sa_mapq == 0: continue
                    
                    if is_te(read.reference_name) != is_te(sa_contig):
                        pos = self.calculate_split_position(read, sa_entry)
                        if pos is None: continue
                        
                        if is_te(read.reference_name):
                            # TE→non-TE split: check against SA contig's length
                            ref_length = self.bam.get_reference_length(sa_contig)
                        else:
                            # non-TE→TE split: check against current reference
                            ref_length = self.bam.get_reference_length(read.reference_name)

                        if pos < 0 or pos >= ref_length:
                            continue
                        
                        # Proceed with valid pos
                        if is_te(read.reference_name):
                            chrom = sa_contig
                            te_id = get_te_family(read.reference_name)
                            # Get correct reference length for non-TE chromosome
                            ref_length = self.bam.get_reference_length(sa_contig)
                        else:
                            chrom = read.reference_name
                            te_id = get_te_family(sa_contig)
                            ref_length = self.bam.get_reference_length(read.reference_name)

                        # Validate position against correct reference
                        if pos < 0 or pos >= ref_length:
                            continue
                        
                        il_candidates[chrom][te_id].append((max(0, pos - 20), min(pos + 20, ref_length)))

        # Merge and process candidates
        for chrom in il_candidates:
            for te_id in il_candidates[chrom]:
                merged = self.merge_intervals(il_candidates[chrom][te_id])
                for start, end in merged:
                    is_novel = True
                    if self.bed_regions:
                        for bed_start, bed_end, bed_family, bed_type in self.bed_regions.get(chrom, []):
                            if te_id == bed_family and not (end < bed_start or start > bed_end):
                                is_novel = False
                                break
                    if is_novel:
                        self.merged_ils[chrom].append((start, end, te_id, 'Novel'))

        # Add BED entries
        if self.bed_regions:
            for chrom in self.bed_regions:
                for start, end, te_family, bed_type in self.bed_regions[chrom]:
                    self.merged_ils[chrom].append((start, end, te_family, bed_type))

    def merge_intervals(self, intervals):
        sorted_intervals = sorted(intervals, key=lambda x: x[0])
        merged = []
        for interval in sorted_intervals:
            if not merged or interval[0] > merged[-1][1] + 20:
                merged.append(list(interval))
            else:
                merged[-1][1] = max(merged[-1][1], interval[1])
        return [tuple(i) for i in merged]

    def save_results(self, output_file):
        with open(output_file, 'w') as f:
            for chrom in self.merged_ils:
                for start, end, te_id, il_type in self.merged_ils[chrom]:
                    f.write(f"{chrom}\t{start}\t{end}\t{te_id}\t{il_type}\n")

if __name__ == '__main__':
    args = parse_args()
    processor = FirstPassProcessor(args.bam, args.bed)
    processor.process_reads()
    processor.save_results(args.output)