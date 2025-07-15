####################################################################
# Microbiome Analysis in R: Compositional Profiles
# Author: Lucky
# Date: [2025-07]
# Description:
#   This script extracts ASV count data from a phyloseq object,
#   summarizes taxonomic composition at multiple levels 
#   (Phylum to Genus), merges with sample metadata, and 
#   computes both raw counts and relative abundances.
#   Output tables are exported to Excel files for further analysis.
####################################################################


# Load workspace and packages
setwd("~/Documents/Obese_Microbiome/4.Compositional_Profiles")
load("~/Documents/Obese_Microbiome/1.Raw_Data/Building_Phyloseq.RData")

library(phyloseq)
library(dplyr)
library(ggplot2)
library(car)  # For Export()

# Helper function: summarize and export by taxonomic level
process_taxonomic_level <- function(level, prefix) {
  message("\nProcessing: ", level)

  # Assign taxonomy
  asv_count_table[[level]] <- taxonomy_table[[level]][match(rownames(asv_count_table), rownames(taxonomy_table))]

  # Group by taxonomy level and sum
  tax_table_level <- asv_count_table %>%
    group_by(.data[[level]]) %>%
    summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(TotalCounts = rowSums(across(where(is.numeric)))) %>%
    arrange(desc(TotalCounts)) %>%
    select(-TotalCounts)

  # Transpose for merging
  tax_table_t <- as.data.frame(t(tax_table_level[,-1]))
  colnames(tax_table_t) <- tax_table_level[[level]]
  tax_table_t$SampleID <- rownames(tax_table_t)

  # Handle NA column names
  colnames(tax_table_t)[is.na(colnames(tax_table_t))] <- "unclassified"

  # Merge with metadata
  tax_with_metadata <- sample_metadata %>%
    left_join(tax_table_t, by = "SampleID")

  # Remove prefix (e.g., "p__")
  prefix_pattern <- paste0("^", substr(level, 1, 1), "__")
  colnames(tax_with_metadata) <- gsub(prefix_pattern, "", colnames(tax_with_metadata))

  # Export absolute counts
  Export(tax_with_metadata, paste0(prefix, "_Counts_", level, ".xlsx"), row.names = FALSE)

  # Compute relative abundance
  tax_with_metadata_rel <- tax_with_metadata %>%
    rowwise() %>%
    mutate(across(
      .cols = -(SampleID:Group),
      .fns = ~ .x / sum(c_across(-(SampleID:Group))),
      .names = "{.col}"
    )) %>%
    ungroup()

  # Export relative abundances
  Export(tax_with_metadata_rel, paste0(prefix, "_Rel_", level, ".xlsx"), row.names = FALSE)
}

# Extract base tables
asv_count_table <- as.data.frame(otu_table(ps))
taxonomy_table <- as.data.frame(tax_table(ps))
sample_metadata <- as(sample_data(ps), "data.frame")

# Sort ASVs by total count
asv_count_table$TotalCounts <- rowSums(asv_count_table)
asv_taxonomy_table_sorted <- asv_count_table[order(-asv_count_table$TotalCounts), ]
rownames(asv_taxonomy_table_sorted) <- NULL
Export(asv_taxonomy_table_sorted, "asv_taxonomy_table_sorted.xlsx")

# Process each taxonomic rank
levels <- c("Phylum", "Class", "Order", "Family", "Genus")
prefixes <- c("1", "2", "3", "4", "5")

for (i in seq_along(levels)) {
  process_taxonomic_level(levels[i], prefixes[i])
}
