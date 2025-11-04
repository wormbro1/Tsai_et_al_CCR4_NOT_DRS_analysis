#!/usr/bin/python
# This script averages the polyA tail length of all reads mapped to each gene
# and includes both valid polyA read count and total read count per gene.

import sys
import pandas as pd

# Read input file
if len(sys.argv) == 2:
    featureCounts_file = sys.argv[1]
else:
    print("Usage: input a featureCounts file that has a polyA tail lengths column")
    sys.exit(1)

# Load file
df = pd.read_csv(featureCounts_file, sep="\t", low_memory=False, na_values=[])

# Ensure polyA_length is numeric; invalid values become NaN
df['polyA_length'] = pd.to_numeric(df['polyA_length'], errors='coerce')

# Group by gene and calculate:
# - mean polyA length (ignores NaN)
# - count of non-NaN polyA lengths
# - total count of reads per gene (including NaNs)
gene_polyA_stats = df.groupby('targets').agg(
    mean_polyA_length=('polyA_length', 'mean'),
    reads_with_polyA=('polyA_length', 'count'),
    total_reads=('polyA_length', 'size')
).reset_index()

# Optional: fill NaN in mean_polyA_length with string 'NA' for clarity
gene_polyA_stats['mean_polyA_length'] = gene_polyA_stats['mean_polyA_length'].round(2)
gene_polyA_stats['mean_polyA_length'] = gene_polyA_stats['mean_polyA_length'].fillna('NA')

# Write result
output_file = featureCounts_file + "_polyAmean.featureCounts"
gene_polyA_stats.to_csv(output_file, sep="\t", index=False)

print(f"All done. Output written to {output_file}")

