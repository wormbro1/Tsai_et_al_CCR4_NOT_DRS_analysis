#!/usr/bin/env bash
#BSUB -J nnpr_pipel_rep1
#BSUB -n 20
#BSUB -q rna
#BSUB -R "select[mem>50] rusage[mem=50] span[hosts=1]"
#BSUB -o logs/z_std_%J.out
#BSUB -e logs/z_std_%J.err

source /beevol/home/ouyangjo/.bashrc
source activate test
module load samtools

##inputing the dorado output bam file path and gtf path (and python script):
bam_file="bam_file_input"
gtf_file="/beevol/home/ouyangjo/genome_files/hg38.ncbiRefSeq_rep1.gtf"
featureCounts_polyAmean_py_script="/beevol/home/ouyangjo/featureCounts_polyAmean_rep1.py"
gene_standardization_script="/beevol/home/ouyangjo/genome_files/gene_name_standardization_rep1.sh"
gene_names="/beevol/home/ouyangjo/genome_files/hgnc_complete_set_rep1.txt"
RPM_calculator="/beevol/home/ouyangjo/RPM_calculator_rep1.sh"

##creating the file name variables that will be used below:
bam_primary_file="${bam_file%.bam}_primary.bam"
sam_file="${bam_file%.bam}_primary_pA.sam"
sam_file_nochrM="${bam_file%.bam}_primary_pA_nochrM.sam"
ENO2_sam="${bam_file%.bam}_primary_pA_ENO2.sam"
ENO2_pA_file="${bam_file%.bam}_primary_pA_ENO2.txt"
featureCounts_gene_assignments_file="${bam_file%.bam}_primary_pA_nochrM_gene_assignments.txt" 
featureCounts_file="${bam_file%.bam}_primary_pA_nochrM.sam.featureCounts"
featureCounts_pAmean_file="${bam_file%.bam}_primary_pA_nochrM.sam.featureCounts_polyAmean.featureCounts"
featureCounts_pAmean_file_RPM="${bam_file%.bam}_primary_pA_nochrM.sam.featureCounts_polyAmean.featureCounts_withRPM.txt"

##filtering the bam file for only primary alignments
samtools view -bF 0x900 -q 1 "$bam_file" > "$bam_primary_file"

##changing the alignment file to a sam file
samtools view -h --threads 20 -o "$sam_file" \
 "$bam_primary_file"
echo "Finished making bam file a sam file!"

##appending the polyA column (pt:t:) to the read names of the sam file:
awk 'BEGIN { OFS = "\t" }
{
    pt_col = "";
    for (i = 2; i <= NF; i++) {
        if ($i ~ /^pt/) {
            pt_col = $i;
            break;
        }
    }
    if (pt_col != "") {
        $1 = $1 "_" pt_col;
    }
    print $0;
}' "$sam_file" \
> tmp && mv tmp "$sam_file"
echo "Finished appending the polyA data to the end of the read name!"

##extracting the ENO2 reads
awk '$3 == "yeast_ENO2"' "$sam_file" > "$ENO2_sam"

##creating the ENO2 file with the polyA tail lengths per read
awk 'BEGIN {
    OFS = "\t"
    print "read_ID", "polyA_length"
}
$2 == 0 {
    read_id = $1
    if (match(read_id, /pt:i:([0-9]+)/, m)) {
       	polyA = m[1]
    } else {
	polyA = "NA"
    }
    print read_id, polyA
}' "$ENO2_sam" > "$ENO2_pA_file"

##removing reads that map to chrM and yeast_ENO2
awk 'BEGIN {OFS="\t"} /^@/ || $3 != "chrM"' "$sam_file" > tmp && mv tmp "$sam_file_nochrM"
awk 'BEGIN {OFS="\t"} /^@/ || $3 != "yeast_ENO2"' "$sam_file_nochrM" > tmp && mv tmp "$sam_file_nochrM"

##do feature counts
featureCounts \
  -a "$gtf_file" \
  -o "$featureCounts_gene_assignments_file"  \
  --largestOverlap \
  -L \
  -R CORE \
  "$sam_file_nochrM"
echo "Finished featureCounts!"

##modify feature counts to make a polyA column
awk 'BEGIN { OFS = "\t" }
{
    match($1, /pt:i:([0-9]+)/, m)
    if (m[1] != "") {
        print $0, m[1]
    } else {
        print $0, "NA"
    }
}' "$featureCounts_file" > tmp && mv tmp "$featureCounts_file"
awk 'BEGIN { print "read_name\tstatus\tnumb_of_targets\ttargets\tpolyA_length" } { print }' \
"$featureCounts_file" > tmp && mv tmp "$featureCounts_file"
echo "Finished creating the polyA_length column in the featureCounts file!"

##condense the feature counts with custom python script
python "$featureCounts_polyAmean_py_script" "$featureCounts_file"
echo "Finished calculating the average reads per gene!"

##calculating RPM
bash "$RPM_calculator" "$featureCounts_pAmean_file" "$featureCounts_file"

##standardizing the gene names:
bash "$gene_standardization_script" "$featureCounts_file" "$gene_names"
bash "$gene_standardization_script" "$featureCounts_pAmean_file" "$gene_names"
bash "$gene_standardization_script" "$featureCounts_pAmean_file_RPM" "$gene_names"

echo "all done!"
