############################################################
# Microbiome Analysis in R: Building a Phyloseq Object
# Author: Lucky
# Date: [2025-07]
# Description:
#   This script loads ASV count data, taxonomy, metadata,
#   and a phylogenetic tree to build a Phyloseq object.
#   For QIIME2 outputs processed for use in R.
############################################################


# ------------------------------
# 1. Install & Load Packages
# ------------------------------

# Install Bioconductor manager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


# Install Bioconductor packages
BiocManager::install("phyloseq")
BiocManager::install("biomformat") # to read and manipulate BIOM files.


# Install CRAN packages
install.packages(c("ape", "dplyr", "readxl", "tibble", "rio", "car"))


# Load packages
library(phyloseq)
library(ape)
library(readxl)
library(dplyr)
library(tibble)
library(rio)
library(car)




# ------------------------------
# 2. Import Input Files
# ------------------------------

# Set working directory (you can change your directory name) 
setwd("~/Documents/Microbiome_Analysis_R/1.Raw_Data")

# Define input file paths
asv_table_file <- "feature-table.xlsx"    # ASV/OTU count table
taxonomy_file  <- "taxonomy.tsv"          # Taxonomy table
metadata_file  <- "Metadata.txt"          # Sample metadata
tree_file      <- "tree.nwk"              # Rooted phylogenetic tree

# Load data
feature_table   <- read_excel(asv_table_file)
tax_table_raw   <- read.delim(taxonomy_file, sep = "\t", header = TRUE)
metadata        <- read.delim(metadata_file)
phy_tree        <- read_tree(tree_file)


# Exploring data
head(feature_table)
head(tax_table_raw)
head(metadata)




# ------------------------------
# 3. Prepare Feature Table
# ------------------------------

## Data manipulation before constructing the phyloseq object

# Phyloseq objects need to have row.names. Define the row names from the ASV_ID column
feature_table_mat <- feature_table %>%
  column_to_rownames("ASV_ID")

head(feature_table_mat)


# Check sample order in metadata to match the sample ID in feature_table_mat
sample_order <- names(feature_table_mat)
sample_order

# To check if SampleID are in the same order
identical(sample_order, metadata$SampleID)
# [FALSE]


# [FALSE] then, Reorders metadata to have rows in the same order as sample_order, based on matching the sampleid column.
metadata_matched <- metadata %>%
  slice(match(sample_order, sampleid)) %>% 
  rename(SampleID = sampleid,
         Group = Group) 

head(metadata_matched)


# Set rownames as SampleID
rownames(metadata_matched) <- metadata_matched$SampleID
head(metadata_matched)


# To check if SampleID are in the same order
identical(sample_order, metadata_matched$SampleID)
# [TRUE] You can proceed the next step.




# ------------------------------
# 4. Prepare Taxonomy Table
# ------------------------------

# Format taxonomy table for phyloseq
# Keep only necessary columns
taxonomy <- tax_table_raw[, c("Feature.ID", "Taxon")]
taxonomy$Taxon <- as.character(taxonomy$Taxon)
head(taxonomy)


# Split taxonomy string into ranks
taxonomy_split <- strsplit(taxonomy$Taxon, "; ")
head(taxonomy_split)


# Create taxonomy matrix (6 levels: Kingdom to Genus)
taxonomy_matrix <- do.call(rbind, lapply(taxonomy_split, function(x) {
  length(x) <- 6
  return(x)
}))

# Assign row names as ASV IDs
rownames(taxonomy_matrix) <- taxonomy$Feature.ID
head(taxonomy_matrix)




# ------------------------------
# 5. Build the Phyloseq Object
# ------------------------------

# Build the phyloseq object
ps <- phyloseq(otu_table(feature_table_mat,  taxa_are_rows = TRUE), 
               tax_table(taxonomy_matrix), 
               sample_data(metadata_matched), 
               phy_tree(phy_tree))

# Confirm structure
ps



# Rename taxonomic ranks inside your tax_table
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", 
                             "Order", "Family", "Genus")
rank_names(ps)




# -----------------------------------------
# 6. Rarefying for determining diversity
# ----------------------------------------

# Check for rarefaction:
library_sizes <- sample_sums(ps)
summary(library_sizes)

# Rarefying
ps_rare <- rarefy_even_depth(
  ps,
  sample.size = min(sample_sums(ps)),  # min sequencing depth
  rngseed = 123,                       # correct argument for reproducibility
  replace = FALSE,                    # sampling without replacement
  verbose = TRUE
)

# Confirm that all samples now have the same depth:
summary(sample_sums(ps_rare))




# ------------------------------
# 7. Explore Phyloseq Object
# ------------------------------

# Basic inspection
sample_names(ps_rare)[1:5]
tax_table(ps_rare)[1:5, ]
otu_table(ps_rare)[1:5, ]
rank_names(ps_rare)
sample_variables(ps_rare)
summary(sample_sums(ps_rare))

# Check for missing taxonomy
table(is.na(tax_table(ps_rare)))



# ------------------------------
# 8. Summarize Taxonomy Table
# ------------------------------

# Convert taxonomy table to data frame
taxa_df <- as.data.frame(tax_table(ps))

# Count unique taxa at each rank
summary_table <- data.frame(
  Rank = colnames(taxa_df),
  Unique_Taxa = sapply(taxa_df, function(x) length(unique(na.omit(x))))
)
  
# View summary
print(summary_table)


  

# ------------------------------
# 9. Export Taxonomy Tables
# ------------------------------

# Add ASV ID to taxonomy
taxa_df_asv <- data.frame(ASV_ID = rownames(taxa_df), taxa_df)

# Combine with abundance data
taxa_rank_asv <- data.frame(
  ASV_ID = rownames(taxa_df),
  taxa_df,
  as.data.frame(otu_table(ps))
)

# Export as Excel files;  You can rename the files by changing the string ("taxa_df_asv.xlsx" → "your_custom_filename.xlsx")
Export(taxa_df_asv, "taxa_df_asv.xlsx")
Export(taxa_rank_asv, "taxa_rank_asv.xlsx")
