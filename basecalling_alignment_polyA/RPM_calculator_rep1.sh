#!/bin/bash

# Usage: ./add_RPM.sh file1 file2
# file1: input featureCounts-like file (tab-delimited, must contain "total_reads" column)
# file2: reference file (used to count rows for normalization)

if [ $# -ne 2 ]; then
  echo "Usage: $0 file1 file2"
  exit 1
fi

file1=$1
file2=$2

# Number of rows in file2 minus 1 for header
rows=$(($(wc -l < "$file2") - 1))

# Get column index of "total_reads"
col=$(head -1 "$file1" | tr '\t' '\n' | grep -nx "total_reads" | cut -d: -f1)

if [ -z "$col" ]; then
  echo "Error: total_reads column not found in $file1"
  exit 1
fi

# Add RPM column: (total_reads / rows) * 1,000,000
awk -v col="$col" -v rows="$rows" 'BEGIN{OFS="\t"}
NR==1 {print $0,"RPM"; next}
{
  rpm = ($col / rows) * 1000000
  print $0, rpm
}' "$file1" > "${file1%.txt}_withRPM.txt"

echo "Done. Output written to ${file1%.txt}_withRPM.txt"

