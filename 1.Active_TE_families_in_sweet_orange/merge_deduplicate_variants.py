import sys

# A function to load variants from a file
def load_variants(filename):
    variants = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('V'):  # Ensure it's a variant line
                parts = line.strip().split()
                # Convert relevant fields to integers for numerical operations
                chrom, start, end, qstart, qend = parts[1], int(parts[2]), int(parts[3]), int(parts[9]), int(parts[10])
                ref_len = end - start
                query_len = qend - qstart
                variant_type = 'INS' if query_len > 50 and ref_len < 20 else 'DEL' if query_len < 20 and ref_len >= 50 else 'OTHER'
                variants.append((chrom, start, end, variant_type, query_len, line))
    return variants

# A function to check if two variants should be considered the same
def is_same_variant(v1, v2):
    # Check if within 10 bp sway in indel end coordinates
    if abs(v1[1] - v2[1]) <= 10 and abs(v1[2] - v2[2]) <= 10:
        # Insertion criteria
        if v1[3] == 'INS' and v2[3] == 'INS':
            return abs(v1[4] - v2[4]) / max(v1[4], v2[4]) < 0.05
        # Deletion criteria
        elif v1[3] == 'DEL' and v2[3] == 'DEL':
            return True
    return False

# Merge and deduplicate variants
def merge_variants(variants):
    merged = []
    for v in variants:
        if not merged:
            merged.append(v)
        else:
            if not any(is_same_variant(v, mv) for mv in merged):
                merged.append(v)
    return merged

# Main script starts here
filenames = sys.argv[1:]  # Get .var.txt filenames from command line arguments
all_variants = []
for filename in filenames:
    all_variants.extend(load_variants(filename))

all_variants.sort(key=lambda x: (x[0], x[1]))  # Sort by chromosome and start
deduplicated_variants = merge_variants(all_variants)

# Output deduplicated variants
for v in deduplicated_variants:
    print(v[5].strip())  # v[5] contains the original line from the file
