############################################################
# Microbiome Analysis in R: Building a Phyloseq Object
# Author: Lucsame Gruneck (Lucky)
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
install.packages(c("ape", "dplyr", "readxl", "tibble"))


# Load packages
library(phyloseq)
library(ape)
library(readxl)
library(dplyr)
library(tibble)



# ------------------------------
# 2. Import Input Files
# ------------------------------

# Set working directory
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
phy_tree_object <- read_tree(tree_file)


# Exploring data
head(feature_table)
head(tax_table_raw)
head(metadata)




#---------------------------------
# 3. Building the Phyloseq Object
#---------------------------------

## Data manipulation before constructing the phyloseq object

# Phyloseq objects need to have row.names. Define the row names from the ASV_ID column
feature_table_mat <- feature_table %>%
  tibble::column_to_rownames("ASV_ID") 

head(feature_table_mat)

# Check sample order in metadata to match the sample ID in feature_table_mat
sample_order <- names(feature_table_mat)
sample_order

# To check if SampleID are in the same order
identical(sample_order, metadata$SampleID)
# [FALSE]

# Reorders metadata to have rows in the same order as sample_order, based on matching the sampleid column.
metadata_matched <- metadata %>%
  slice(match(sample_order, sampleid)) %>% 
  rename(SampleID = sampleid,
         Group = Group) 

head(metadata_matched)

# Add the SampleID column to rows of metadata_matched
rownames(metadata_matched) <- metadata_matched$SampleID

head(metadata_matched)

# To check if SampleID are in the same order
identical(sample_order, metadata_matched$SampleID)
# [TRUE]


# Format taxonomy table for phyloseq
taxonomy <- tax_table_raw
head(taxonomy)
taxonomy <- taxonomy[, c("Feature.ID", "Taxon")] #Select only "Feature.ID" and "Taxon" columns
taxonomy$Taxon <- as.character(taxonomy$Taxon)
taxonomy_split <- strsplit(taxonomy$Taxon, "; ") # Split taxonomic ranks

head(taxonomy_split)

# Create a matrix with columns: Kingdom, Phylum, Class, Order, Family, Genus
taxonomy_matrix <- do.call(rbind, lapply(taxonomy_split, function(x) {
  length(x) <- 6  # Ensure 6 levels
  return(x)
}))

# Add the Feature.ID column to rows of taxonomy_matrix
rownames(taxonomy_matrix) <- taxonomy$Feature.ID

head(taxonomy_matrix)


# Build the phyloseq object
ps <- phyloseq(otu_table(feature_table_mat,  taxa_are_rows = TRUE), 
               tax_table(taxonomy_matrix), 
               sample_data(metadata_matched), 
               phy_tree(phy_tree))

ps

# Exploring Phyloseq Components
sample_names(ps)[1:10]
tax_table(ps)[1:10]
otu_table(ps)[1:10]
rank_names(ps)
sample_variables(ps)
summary(sample_sums(ps))

#  rename the ranks inside your tax_table
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", 
                             "Order", "Family", "Genus")
rank_names(ps)


# If you want to check if there are any missing classifications, you can do:
table(is.na(tax_table(ps)))


# Get OTU counts at each taxonomy level
## Make a data frame showing each ASV's taxonomy
taxa_df <- as.data.frame(tax_table(ps))

## Summarize the number of unique ASVs at each taxonomic rank (exclude unclassified)
summary_table <- data.frame(
  Rank = c("Kingdon", "Phylum", "Class", "Order", "Family", "Genus"),
  Unique_Taxa = sapply(taxa_df, function(x) length(unique(na.omit(x))))
)

summary_table

# add ASV column to taxa_df
taxa_df_asv <- data.frame(ASV_ID = rownames(taxa_df),
                          taxa_df)
taxa_rank_asv <- data.frame(ASV_ID = rownames(taxa_df),
                            taxa_df, data.frame(otu_table(ps))
                            )

Export(taxa_df_asv, "taxa_df_asv.xlsx")
Export(taxa_rank_asv, "taxa_rank_asv.xlsx")
