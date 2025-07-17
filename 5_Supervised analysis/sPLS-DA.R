#################################################################################
# Microbiome Analysis in R: Metabolite Classification (Positive & Negative Ions)
# Author: Lucky
# Date: 2025-7
# Description:
#   1. Load and prepare positive/negative ion metabolite datasets
#   2. Compute log2 fold change (log2FC) relative to the Normal group
#   3. Summarize group-wise mean ± SD of log2FC for each metabolite
#   4. Visualize results using:
#        - Heatmaps of log2FC and log2-transformed relative abundance
#        - Grouped bar plots colored by direction of regulation
#   5. Export summary statistics for downstream analysis
#   6. Perform statistical analysis (direct comparison)
################################################################################


# Set Working Directory
setwd("~/Documents/Microbiome_Analysis_R/5.sPLS-DA")

# Load Required Libraries
library(readxl)
library(mixOmics)
library(MASS)
library(lattice)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Load Data
load("~/Documents/Microbiome_Analysis_R/5.Metabolites/Metabolites.RData")
metabo_descp <- read_excel("~/Documents/Microbiome_Analysis_R/4.Metabolites/metabo_descp.xlsx")









