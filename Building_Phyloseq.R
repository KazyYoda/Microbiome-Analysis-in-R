################################################
# "Microbiome Analysis: Building Phyloseq in R #
################################################


#--------------------------------
# 1. Loading Required Packages
#--------------------------------

## Install Bioconductor manager if it's not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

## Install biomformat from Bioconductor
BiocManager::install("biomformat") # to read and manipulate BIOM files.

## Install other packages:
install.packages("ape") # Analyses of Phylogenetics and Evolution
install.packages("dplyr") # Data manipulation
install.packages("readxl") # Read Excel Files




#--------------------------------
# 2. Importing Input Files
#--------------------------------

# import input files:
# •	ASV/OTU Table
# •	Taxonomy Table
# •	Sample Metadata
# •	Phylogenetic Tree


## Prepare each datasetset for use in building the phyloseq object
# Set a working directory:
setwd("~/Documents/Microbiome_Analysis_R/1.Raw_Data")


# Load packages:
library(phyloseq)
library(ape)
library(readxl)
library(dplyr)

# New Metadata
BMI_zscore <- final_bmi_for_age %>%
  select(Sample, final_z.indv, BMI_for_age) %>%
  rename("BMI_zscore" = "final_z.indv",
         "Group" = "BMI_for_age")

BMI_zscore$Group[BMI_zscore$Group == "Thinness"] <- "N"
BMI_zscore$Group[BMI_zscore$Group == "Normal"] <- "N"
BMI_zscore$Group[BMI_zscore$Group == "OV"] <- "OW"


# Relabel incorrect BMI groups for sample id BH299 & BH325 in metadata
metadata$Group1[metadata$sampleid == "BH299"] <- "N"
metadata$Group1[metadata$sampleid == "BH325"] <- "OB"

# Match BMI_zscore to Metadata
Metadata_BMIzscore <- metadata %>%
  left_join(BMI_zscore, by = c("sampleid" = "Sample",
                               "Group1" = "Group")) %>%
  rename("Group" = "Group1") %>%
  select(-Group2)


Metadata_new <- Metadata_BMIzscore %>%
  select(sampleid, Group)

Export(Metadata_new, "Metadata_new.txt")


# Open feature-table.tsv with Excel and remove the first row "Constructed from biom file"
# Then, save it as .xlsx file and import in R
# File paths
asv_table_file <- "feature-table.xlsx"
taxonomy_file <- "taxonomy.tsv"
metadata_file <- Metadata_new
tree_file <- "tree.nwk"

# Load data
feature_table <- read_excel("feature-table.xlsx")
tax_table_raw <- read.delim("taxonomy.tsv", sep = "\t", header = TRUE)
metadata <- metadata_file
phy_tree <- read_tree(tree_file)  # if using phylogenetic tree


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



