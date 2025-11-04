# Load libraries
library(DESeq2)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)

rm(list=ls())
directory <- "/Users/jouyang3/Documents/A_Rissland_lab/Sequence_analysis/JiYoun_data/all_reps/DESeq2_analysis/"
setwd(directory)

# 1. Read in sample info (file_map.csv)
file_map <- read.csv("file_map.csv", stringsAsFactors = FALSE)

# 2. Build list of count data frames
counts_list <- lapply(1:nrow(file_map), function(i) {
  f <- file_map$file[i]
  sample_name <- file_map$sample[i]
  
  # Read file safely (allow extra columns)
  df <- read.delim(
    f,
    sep = "\t",
    header = TRUE,
    fill = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Inspect available columns if needed
  # print(colnames(df))
  
  # Select only the gene and count columns
  df <- df[, c("standardized_gene_name_noNA", "total_reads")]
  
  # Rename to "gene" and sample name
  colnames(df) <- c("gene", sample_name)
  return(df)
})

# 3. Merge all count data frames by gene
counts_merged <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), counts_list)
counts_merged[is.na(counts_merged)] <- 0   # replace NAs with 0

# 4. Build DESeq2 count matrix
count_matrix <- as.matrix(counts_merged[,-1])
rownames(count_matrix) <- counts_merged$gene

# 5. Build sample metadata
# Assumes sample names contain condition in their name (control, asv, asv_lipo)
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(gsub("_rep[0-9]+", "", colnames(count_matrix))) # e.g., "control", "asv", "asv_lipo"
)
# Set "control" as the baseline condition
sample_info$condition <- relevel(sample_info$condition, ref = "control")

# 6. Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)
dds <- DESeq(dds)

# -----------------------
# 4. Normalized counts
# -----------------------
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts) %>%
  rownames_to_column("gene")

# -----------------------
# 5. Function to build contrast-specific table
# -----------------------
make_contrast_df <- function(dds_obj, norm_counts_df, sample_info, cond1, cond2) {
  # Get DESeq2 results
  res <- results(dds_obj, contrast = c("condition", cond1, cond2)) %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  # Get sample names for each condition
  samples_cond1 <- rownames(sample_info)[sample_info$condition == cond1]
  samples_cond2 <- rownames(sample_info)[sample_info$condition == cond2]
  
  # Subset normalized counts for these samples
  counts_subset <- norm_counts_df[, c("gene", samples_cond1, samples_cond2)]
  
  # Merge DESeq2 stats with normalized counts
  merged_df <- res %>%
    left_join(counts_subset, by = "gene")
  
  # Add average RPM per condition
  merged_df <- merged_df %>%
    rowwise() %>%
    mutate(
      avg_RPM_cond1 = mean(c_across(all_of(samples_cond1)), na.rm = TRUE),
      avg_RPM_cond2 = mean(c_across(all_of(samples_cond2)), na.rm = TRUE)
    ) %>%
    ungroup()
  
  return(merged_df)
}


# -----------------------
# 6. Create each contrast dataframe
# -----------------------
res_ctrl_asv_df <- make_contrast_df(dds, norm_counts_df, sample_info, "asv", "control")
res_ctrl_asvlipo_df <- make_contrast_df(dds, norm_counts_df, sample_info, "asv_lipo", "control")
res_asv_asvlipo_df <- make_contrast_df(dds, norm_counts_df, sample_info, "asv_lipo", "asv")

# -----------------------
# 7. Save results
# -----------------------
write.csv(res_ctrl_asv_df, "DESeq2_asv_vs_control_with_RPM.csv", row.names = FALSE)
write.csv(res_ctrl_asvlipo_df, "DESeq2_asvlipo_vs_control_with_RPM.csv", row.names = FALSE)
write.csv(res_asv_asvlipo_df, "DESeq2_asvlipo_vs_asv_with_RPM.csv", row.names = FALSE)


##making some graphs:
##making some graphs:
# Prepare data (remove NA padj values)
volcano_df <- res_ctrl_asv_df %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  ) %>%
  filter(!is.na(padj))

# Make volcano plot
ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.7, size = 3) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "gray70"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: asv vs control",
    x = "log2(Fold Change)",
    y = expression(-log[10](adjusted~p~value)),
    color = "Category"
  ) +
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )






##For MA plot:
# Filter out rows with NA or zero baseMean
ma_data <- res_ctrl_asv_df %>%
  filter(!is.na(baseMean) & baseMean > 0)

# Make the MA plot
ggplot(ma_data, aes(x = log10(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5, size = 3) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 2) +
  geom_hline(yintercept = c(1, -1), color = "blue", linetype = "dashed", size = 2) +
  labs(
    x = expression(log[10]~"Mean expression (baseMean)"),
    y = expression(log[2]~"Fold Change"),
    color = "Significant (padj < 0.05)",
    title = "MA Plot: asv vs control"
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )



