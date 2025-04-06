#!/usr/bin/env python3

import pysam
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Second pass IL counting')
    parser.add_argument('-b', '--bam', required=True, help='Input BAM file')
    parser.add_argument('-m', '--merged', required=True, help='Merged ILs BED file')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    return parser.parse_args()

def is_te(contig):
    return contig.startswith('TE')

def get_te_family(te_id):
    return te_id.split('.')[0] if '.' in te_id else te_id

class SecondPassProcessor:
    def __init__(self, bam_path, merged_bed):
        self.bam = pysam.AlignmentFile(bam_path, 'rb')
        self.counts = defaultdict(lambda: {'A':0, 'B':0, 'C':0})
        self.il_index = self.load_merged_ils(merged_bed)

    def load_merged_ils(self, bed_file):
        index = defaultdict(lambda: defaultdict(list))
        with open(bed_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                chrom, start, end, te_id, il_type = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
                index[chrom][te_id].append((start, end, il_type))
        return index

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

    def classify_pair(self, reads):
        """Classify read pairs against existing ILs only"""
        # Track counted ILs to avoid double-counting
        counted_ils = set()

        # Type A classification
        for read in reads:
            if read.has_tag('SA'):
                for sa_entry in read.get_tag('SA').split(';'):
                    if not sa_entry:
                        continue

                    # Parse SA fields
                    sa_fields = sa_entry.split(',')
                    if len(sa_fields) < 6:
                        continue

                    # Extract SA components
                    sa_contig = sa_fields[0]
                    sa_pos = int(sa_fields[1])
                    sa_strand = sa_fields[2]
                    sa_cigar = sa_fields[3]
                    sa_mapq = int(sa_fields[4])

                    # Skip low quality
                    if sa_mapq == 0 or read.mapping_quality == 0:
                        continue

                    # Determine TE context
                    ref_is_te = is_te(read.reference_name)
                    sa_is_te = is_te(sa_contig)
                    if ref_is_te == sa_is_te:
                        continue

                    # Calculate split position
                    split_pos = self.calculate_split_position(read, sa_entry)
                    if not split_pos:
                        continue

                    # Determine IL coordinates and family
                    if ref_is_te:
                        il_chrom = sa_contig
                        te_family = get_te_family(read.reference_name)
                    else:
                        il_chrom = read.reference_name
                        te_family = get_te_family(sa_contig)

                    # Update unpacking to include il_type
                    for il_start, il_end, il_type in self.il_index.get(il_chrom, {}).get(te_family, []):
                        if il_start <= split_pos <= il_end:
                            key = (il_chrom, il_start, il_end, te_family)
                            if key not in counted_ils:
                                self.counts[key]['A'] += 1
                                counted_ils.add(key)

        # Type B classification
        ref_contigs = {r.reference_name for r in reads if not r.is_unmapped}
        if len(ref_contigs) == 2:
            te_refs = [ctg for ctg in ref_contigs if is_te(ctg)]
            if len(te_refs) == 1:
                te_family = get_te_family(te_refs[0])
                nonte_chrom = [ctg for ctg in ref_contigs if not is_te(ctg)][0]
                
                # Check all non-TE reads in pair
                for read in reads:
                    if read.mapping_quality == 0:
                        continue
                    if read.reference_name == nonte_chrom:
                        # Update unpacking here too
                        for il_start, il_end, il_type in self.il_index.get(nonte_chrom, {}).get(te_family, []):
                            if abs(read.reference_start - (il_start + il_end)/2) < 1000:
                                key = (nonte_chrom, il_start, il_end, te_family)
                                if key not in counted_ils:
                                    self.counts[key]['B'] += 1
                                    counted_ils.add(key)

        # Type C classification
        if len(reads) == 2 and all(not is_te(r.reference_name) for r in reads):
            if reads[0].is_proper_pair and reads[1].is_proper_pair:
                chrom = reads[0].reference_name
                start = min(r.reference_start for r in reads)
                end = max(r.reference_end for r in reads)
                
                for te_family in self.il_index.get(chrom, {}):
                    for il_start, il_end, il_type in self.il_index[chrom][te_family]:
                        if (start < il_start) and (il_end < end):
                            key = (chrom, il_start, il_end, te_family)
                            if key not in counted_ils:
                                self.counts[key]['C'] += 1
                                counted_ils.add(key)

    def process_reads(self):
        # Original second_pass processing logic
        self.bam.reset()
        current_qname = None
        reads = []
        
        for read in self.bam:
            if read.is_secondary or read.is_supplementary or read.mapping_quality == 0:
                continue
            
            if read.qname != current_qname:
                if reads:
                    self.classify_pair(reads)
                current_qname = read.qname
                reads = [read]
            else:
                reads.append(read)
        
        if reads:
            self.classify_pair(reads)

    def generate_report(self):
        with open(self.args.output, 'w') as f:
            f.write("Chrom\tStart\tEnd\tTE_Family\tType\tTypeA\tTypeB\tTypeC\tPASS_FILTER\n")
            for chrom in self.il_index:
                for te_id in self.il_index[chrom]:
                    for start, end, il_type in self.il_index[chrom][te_id]:
                        counts = self.counts.get((chrom, start, end, te_id), {'A':0, 'B':0, 'C':0})
                        hq = 'Yes' if counts['A'] > 0 and counts['A'] + counts['B'] > 1 else 'No'
                        f.write(f"{chrom}\t{start}\t{end}\t{te_id}\t{il_type}\t"
                                f"{counts['A']}\t{counts['B']}\t{counts['C']}\t{hq}\n")

    def run(self, args):
        self.args = args
        self.process_reads()
        self.generate_report()

if __name__ == '__main__':
    args = parse_args()
    processor = SecondPassProcessor(args.bam, args.merged)
    processor.run(args)