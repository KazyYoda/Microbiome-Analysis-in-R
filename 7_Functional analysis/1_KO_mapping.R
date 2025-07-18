################################################################################
# Microbiome Functional Analysis: KOs Mapping on KEGG Database
# Author: Lucky
# Date: 2025-07
# Description:
#   1. Import unstratified KO abundance data from PICRUSt2 output
#   2. Clean KO identifiers by removing the "ko:" prefix
#   3. Extract KO IDs and export as plain text for KEGG mapping
#   4. KEGG pathway mapping to be done via Python script in Terminal
################################################################################

# Set working directory
setwd("~/Documents/Microbiome_Analysis_R/8.Functional_analysis/1_KO_Mapping")

# Load required libraries
library(dplyr)
library(readr)
library(rio)

# Step 1: Import KO data (unstratified)
KO_unstrat <- read_delim(
  "~/Documents/Microbiome_Analysis_R/8.Functional_analysis/pred_metagenome_unstrat.tsv",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

# Step 2: Remove "ko:" prefix from the function column
KO_unstrat$`function` <- sub("^ko:", "", KO_unstrat$`function`)

# Step 3: Extract KO IDs
KO_ID <- data.frame(KO_ID = KO_unstrat$`function`)

# Step 4: Export KO ID list for use in KEGG mapping (external Python tool)
export(KO_ID, "KO_ID.txt")



# Then perform KEGG Database Mapping using python running in the "Terminal"













