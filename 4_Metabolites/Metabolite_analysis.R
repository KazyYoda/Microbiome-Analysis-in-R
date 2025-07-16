###########################################################################
# Microbiome Analysis in R: Metabolite Profiling (Positive & Negative Ions)
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
###########################################################################


# ============================
# 1. Setup
# ============================

# Set working directory
setwd("~/Documents/Microbiome_Analysis_R/4.Metabolites")

# Load required packages
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggh4x)

# Define group colors
group_colors <- c(N = "grey", OW = "#FFA500", OB = "darkred")




# ============================
# 2. Load and Prepare Data
# ============================

# Load metabolite data (clean data)
pos_metabo <- read_excel("POS_Metabo.xlsx")
neg_metabo <- read_excel("NEG_Metabo.xlsx")

# If data contain 0/NA then add 1 to numeric columns only to avoid log2(0)
data <- data %>% mutate(across(where(is.numeric), ~.x + 1))

# Check SampleID consistency
identical(pos_metabo$SampleID, neg_metabo$SampleID)

# Extract metadata (n = 55)
Metadata_55 <- pos_metabo %>% select(SampleID, Group)
Metadata_55$Group <- factor(Metadata_55$Group, levels = c("N", "OW", "OB"))

# Apply group as factor
pos_metabo$Group <- factor(pos_metabo$Group, levels = c("N", "OW", "OB"))
neg_metabo$Group <- factor(neg_metabo$Group, levels = c("N", "OW", "OB"))

## Assign coding to both positive and negative metabolite data and remove metadata columns for downsteam analysis
# Preparing POS metabo data
posmetabo <- pos_metabo %>%
  rename_with(
    .cols = 3:97,
    .fn = ~ paste0("pos", seq_along(.))
  ) %>% 
  select(-SampleID, -Group) %>%
  as.data.frame

# Assign SampleIDs to rownames
rownames(posmetabo) <- Metadata_55$SampleID


# Preparing NEG metabo data
negmetabo <- neg_metabo %>%
  rename_with(
    .cols = 3:40,
    .fn = ~ paste0("neg", seq_along(.))
  ) %>% 
  select(-SampleID, -Group) %>%
  as.data.frame

# Assign SampleIDs to rownames
rownames(negmetabo) <- Metadata_55$SampleID




# ============================
# 3. Helper Functions
# ============================

# Function to compute log2FC relative to normal group (N)
log2FC_relative_to_N <- function(metabo_data, metadata, control = "N") {

  # Compute control means from Normal group
  control_means <- colMeans(metabo_data[metadata$Group == control, ], na.rm = TRUE)

  # Use sweep to divide each column by its corresponding control mean (column-wise division) 
  log2FC <- log2(sweep(metabo_data, 2, control_means, "/"))
  return(log2FC)
}


# Function to summarize log2FC by group
summarize_log2FC <- function(log2FC_df, metadata, tag = "pos") {

  # Combine with Metadata
  combined <- data.frame(metadata, log2FC_df)

  # Reshape to long format
  log2FC_long <- pivot_longer(combined, cols = starts_with(tag), names_to = "Metabolite", values_to = "log2FC")

  # Calculate mean log2FC and SD per group 
  summary <- log2FC_long %>%
    group_by(Group, Metabolite) %>%
    summarise(mean_log2FC = mean(log2FC, na.rm = TRUE),
              sd_log2FC = sd(log2FC, na.rm = TRUE), .groups = "drop") %>%
    filter(Group != "N") %>%

    # Classify direction
    mutate(Direction = case_when(
      mean_log2FC > 0 ~ "Upregulated",
      mean_log2FC < 0 ~ "Downregulated",
      TRUE ~ "No change"
    ))
  return(summary)
}


# Function to create bar plot
plot_log2FC_bar <- function(summary_df, title, xlab, fill_colors) {

  # Count number of up/downregulated metabolites per group
  summary_df <- summary_df %>%
    arrange(Direction, mean_log2FC) %>%
    mutate(Metabolite = factor(Metabolite, levels = unique(Metabolite)))

  # Define group color
  group_colors <- c(N = "grey", OW = "#FFA500", OB = "darkred")


  # Create a bar plot of log2FC with error bars (±SD), colored by direction, and ordered by mean_log2FC
  p <- ggplot(summary_df, aes(x = Metabolite, y = mean_log2FC, fill = Direction)) +
    geom_col(position = position_dodge(0.8), width = 0.7, aes(group = Group)) +
    geom_text(aes(label = round(sd_log2FC, 2),
                  y = ifelse(mean_log2FC >= 0, mean_log2FC + 0.8, mean_log2FC - 0.8)),
              size = 1.5, position = position_dodge(0.8),
              hjust = ifelse(summary_df$mean_log2FC >= 0, 0, 1)) +
    scale_fill_manual(values = fill_colors) +
    facet_wrap2(~ Group, strip = strip_themed(
      background_x = elem_list_rect(fill = group_colors[2:3], color = "black"),
      text_x = elem_list_text(face = "bold", color = c("black", "white"))
    )) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 8, face = "bold"),
          strip.text = element_text(face = "bold")) +
    labs(x = xlab, y = "Mean log2 Fold Change", fill = "Direction", title = title) +
    coord_flip()

  return(p)
}


# Function to create heatmap
plot_heatmap <- function(mat, metadata, legend_title = "Log2") {


  # Ensure column names of the transposed matrix match SampleID 
  colnames(mat) <- metadata$SampleID

  # Print sample dimensions
  cat("Matrix columns (samples):", ncol(mat), "\n")               # e.g., "Matrix columns (samples): 55"
  cat("Metadata samples       :", length(metadata$Group), "\n")   # e.g., "Metadata samples       : 55"

  # Check match
  if (ncol(mat) == length(metadata$Group)) {
    cat("✅ Sample counts matched.\n\n")
  } else {
    cat("❌ Sample counts do NOT match! Please check your data.\n\n")
  }


  # Define group color
  group_colors <- c(N = "grey", OW = "#FFA500", OB = "darkred")

  # Create sample group annotation
  sample_group <- data.frame(Group = Metadata_55$Group)
  rownames(sample_group) <- Metadata_55$SampleID

  # Define the annotation
  ha <- HeatmapAnnotation(df = sample_group, col = list(Group = group_colors))

  # Color function for heatmap
  col_fun <- colorRamp2(c(min(mat), mean(mat), max(mat)), c("#FC766A", "grey96", "#88C0A7"))

  # Create the heatmap with row clustering based on Spearman correlation
  Heatmap(
    mat,
    name = legend_title,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    border = TRUE,
    border_gp = gpar(col = "grey"),
    rect_gp = gpar(col = "grey30", lwd = 0.5),
    col = col_fun,
    top_annotation = ha,
    column_split = metadata$Group,
    column_names_gp = gpar(fontsize = 6),
    row_names_gp = gpar(fontsize = 6),
    clustering_distance_rows = "spearman",
    clustering_distance_columns = "spearman",
    heatmap_legend_param = list(direction = "vertical")
  )
}




# ---------------------------------------------------------------
# Example Usage:
# ---------------------------------------------------------------

# 1. Compute log2FC relative to normal group (N)
log2FC_pos_reltoN <- log2FC_relative_to_N(posmetabo, Metadata_55, control = "N")
log2FC_neg_reltoN <- log2FC_relative_to_N(negmetabo, Metadata_55, control = "N")

Export(log2FC_pos_reltoN, "1_log2FC_pos_reltoN.xlsx")
Export(log2FC_neg_reltoN, "2_log2FC_neg_reltoN.xlsx")


#2.  Summarize log2FC by group
log2FC_grouped_pos <- summarize_log2FC(log2FC_pos_reltoN, Metadata_55, "POS")
log2FC_grouped_neg <- summarize_log2FC(log2FC_neg_reltoN, Metadata_55, "NEG")

Export(log2FC_grouped_pos, "1_log2FC_grouped_pos.xlsx")
Export(log2FC_grouped_neg, "2_log2FC_grouped_neg.xlsx")


# 3. Plot grouped barplot
ggpos <- plot_log2FC_bar(summary_df = log2FC_grouped_pos,
                         title = "Log2FC for Positive-Ion Metabolites (vs. Normal)",
                         xlab = "Metabolite",
                         fill_colors = c("Upregulated" = "#88C0A7", "Downregulated" = "#FC766A"))

ggneg <- plot_log2FC_bar(summary_df = log2FC_grouped_neg,
                         title = "Log2FC for Negative-Ion Metabolites (vs. Normal)",
                         xlab = "Metabolite",
                         fill_colors = c("Upregulated" = "#88C0A7", "Downregulated" = "#FC766A"))

# Print plot
print(ggpos)
print(ggneg)


# 4. Heatmap plot
# Transpose POS log2FC matrix for heatmap (Metabolites as rows)
log_matrix_pos <- as.matrix(log2(posmetabo))
rownames(log_matrix_pos) <- Metadata_55$SampleID
log2_matrix_pos_t <- t(log_matrix_pos)

# Generate heatmap for POS log2FC
plot_heatmap(mat = log2_matrix_pos_t,
             metadata = Metadata_55,
             legend_title = "Log2 (POS)")

plot_heatmap(mat = log2_matrix_neg_t,
             metadata = Metadata_55,
             legend_title = "Log2 (NEG)")
