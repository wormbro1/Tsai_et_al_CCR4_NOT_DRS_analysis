#!/bin/bash

# Usage:
# bash gene_name_standardization.sh input.tsv hgnc_complete_set.txt

file1="$1"
file2="$2"
outfile="${file1%.tsv}_standardized.tsv"

awk -F'\t' -v OFS='\t' -v outfile="$outfile" '
BEGIN {
  IGNORECASE = 1
}

# Parse file2 header to find column indices
FNR == 1 && NR == FNR {
  for (i = 1; i <= NF; i++) {
    col = tolower($i)
    gsub(/"/, "", col)
    if (col == "symbol") sym_col = i
    else if (col == "prev_symbol") prev_col = i
    else if (col == "alias_symbol") alias_col = i
  }
  next
}

# Build alias map from file2
NR == FNR {
  sym = toupper($sym_col)
  gsub(/"/, "", sym)
  if (sym != "") {
    if ((!(sym in priority)) || priority[sym] > 1) {
      alias[sym] = sym
      priority[sym] = 1
    }
  }

  n = split($prev_col, prevs, /\|/)
  for (i = 1; i <= n; i++) {
    prev = toupper(prevs[i])
    gsub(/"/, "", prev)
    if (prev != "" && ((!(prev in priority)) || priority[prev] > 2)) {
      alias[prev] = sym
      priority[prev] = 2
    }
  }

  m = split($alias_col, aliases, /\|/)
  for (i = 1; i <= m; i++) {
    a = toupper(aliases[i])
    gsub(/"/, "", a)
    if (a != "" && ((!(a in priority)) || priority[a] > 3)) {
      alias[a] = sym
      priority[a] = 3
    }
  }
  next
}

# File1 (input): detect target column, apply mapping
NR != FNR {
  if (FNR == 1) {
    for (i = 1; i <= NF; i++) {
      colname = tolower($i)
      gsub(/"/, "", colname)
      if (colname == "targets") target_col = i
      header[i] = $i
    }

    if (!target_col) {
      print "Error: 'targets' column not found in input file" > "/dev/stderr"
      exit 1
    }

    print $0, "standardized_gene_name", "standardized_gene_name_noNA"
    next
  }

  orig = $target_col
  gene = toupper(orig)
  gsub(/"/, "", gene)

  if (gene in alias) {
    std = alias[gene]
  } else {
    std = "NA"
    na++
  }

  std_out = (std == "NA") ? orig : std

  print $0, std, std_out
  total++
}

END {
  print "Total genes processed:\t" total > "/dev/stderr"
  print "Genes without match (NA):\t" na > "/dev/stderr"
  print "Standardized output written to:\t" outfile > "/dev/stderr"

  # Optional: save alias map
  for (a in alias) {
    print a, alias[a] > "_alias_map.tsv"
  }
}
' "$file2" "$file1" > "$outfile"

