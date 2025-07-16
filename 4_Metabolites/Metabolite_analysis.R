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
#   6. Perform statistical analysis
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
data <- data %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0) + 1))

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

# Generate heatmap for POS log2
plot_heatmap(mat = log2_matrix_pos_t,
             metadata = Metadata_55,
             legend_title = "Log2 (POS)")


# Transpose NEG log2FC matrix for heatmap (Metabolites as rows)
log_matrix_neg <- as.matrix(log2(negmetabo))
rownames(log_matrix_neg) <- Metadata_55$SampleID
log2_matrix_neg_t <- t(log_matrix_neg)

# Generate heatmap for NEG log2
plot_heatmap(mat = log2_matrix_neg_t,
             metadata = Metadata_55,
             legend_title = "Log2 (NEG)")




# ===================================================
# 4. Normality and homogeneity of variance checking
# ===================================================

# Function to check normality and homogeneity of variance
norm_var_check <- function(metabo_data, metadata, metabo_descp, tag = "pos") {
    
    # Ensure required library is loaded
    library(dplyr)
    library(car)  # for leveneTest
    
    # Convert input matrix to data frame
    log2_metabo_data <- as.data.frame(metabo_data)
    
    # --- Normality Check using Shapiro-Wilk ---
    shapiro_list <- apply(log2_metabo_data, 2, shapiro.test)
    
    metabo_shapiro <- data.frame(
        Metabolite = metabo_descp$Metabolite,
        Code = metabo_descp$Code,
        W = sapply(shapiro_list, function(x) x$statistic),
        Norm_p.value = sapply(shapiro_list, function(x) x$p.value)
    ) %>%
        mutate(across(where(is.numeric), \(x) round(x, digits = 6)))
    
    # --- Homogeneity of Variance using Levene's Test ---
    levene_list <- lapply(log2_metabo_data, function(x) {
        leveneTest(x ~ factor(metadata$Group), data = metadata)
    })
    
    # Extract p-values from each Levene test
    levene_pvals <- sapply(levene_list, function(result) {
        if (is.data.frame(result) && "Pr(>F)" %in% colnames(result)) {
            return(result$`Pr(>F)`[1])
        } else {
            return(NA)
        }
    })
    
    levene_df <- data.frame(
        Metabolite = names(levene_list),
        Var_p.value = round(levene_pvals, 3)
    )
    
    # --- Combine results ---
    NormVar <- left_join(metabo_shapiro, levene_df, by = "Metabolite")
    return(NormVar)
}

# ---------------------------------------------------------------
# Example Usage: Normality and homogeneity of variance checking
# ---------------------------------------------------------------
pos_norm_var <- norm_var_check(log2_matrix_pos, Metadata_55, pos_descp, "pos")
neg_norm_var <- norm_var_check(log2_matrix_neg, Metadata_55, neg_descp, "pos")


               

# =============================================
# 5. Statistical analysis - direct comparison
# =============================================

# ----- Kruskal-Wallis: Helper function -----
kruskal_metabo <- function(metabo_norm_var, metabo_data, metadata, metabo_descp, tag = "pos") {
    
    # Ensure input matrix is treated as a data frame
    metabo_data <- as.data.frame(metabo_data)
    
    # Step 1: Filter metabolites that fail both normality & equal variance assumptions
    nonNorm <- metabo_norm_var %>% 
        filter(Norm_p.value < 0.05)
    
    # Step 2: Subset non-normal metabolites from the full dataset
    metabo_kruskal <- metabo_data %>% 
        select(all_of(nonNorm$Code))
    
    # Step 3: Perform Kruskal-Wallis test for each metabolite
    krus_results <- lapply(metabo_kruskal, function(x) kruskal.test(x ~ Group, data = metadata))
    
    # Step 4: Extract p-values and format results
    krus_pvalue <- data.frame(
        Code = names(metabo_kruskal),
        p.value = sapply(krus_results, function(res) res$p.value)
    ) %>%
        mutate(p.value = round(p.value, 6))  # Round p-values to 6 digits
    
    # Step 5: Merge with metabolite descriptions
    krus_pvalue_descp <- krus_pvalue %>%
        left_join(nonNorm, by = c("Code" = "Code")) %>%
        select(Code, Metabolite, p.value)
    
    # Print the result to console
    print(krus_pvalue_descp)
}

# ---------------------------------------------------------------
# Example Usage: Direct comparison - kruskal
# ---------------------------------------------------------------
pos_kruskal <- kruskal_metabo(log2_matrix_pos, Metadata_55, pos_descp, "pos")
neg_krukal <- kruskal_metabo(log2_matrix_neg, Metadata_55, neg_descp, "neg")

               

# ----- Dunn's Test - Post-hoc Test for Kruskal-Wallis: Helper function -----
dunnTest_metabo <- function(krus_pvalue_descp, metabo_data, metadata, metabo_descp, tag = "pos") {
    
    # Ensure input matrix is treated as a data frame
    metabo_data <- as.data.frame(metabo_data)
    
    # Rename colnames for metabo_data
    colnames(metabo_data) <- metabo_descp$Metabolite
    
    # Load required package
    library(FSA)
    library(dplyr)
    
    # 1. Filter only metabolites with significant Kruskal-Wallis p-values
    krus_sig <- krus_pvalue_descp %>% 
        filter(p.value < 0.05)
    
    # 2. Run Dunn's test for each significant metabolite
    dunn_test <- setNames(
        lapply(krus_sig$Metabolite, function(x) {
            dunn_res <- dunnTest(metabo_data[[x]] ~ Group, data = metadata, method = "bh")  # BH adjustment
            dunn_df <- as.data.frame(dunn_res$res) 
            return(dunn_df)
        }),
        krus_sig$Metabolite 
    )
    
    # 3. Extract significant pairwise comparisons (adjusted p < 0.05)
    dunn_sig <- do.call(rbind, lapply(names(dunn_test), function(metabolite) {
        df <- dunn_test[[metabolite]]
        df$Metabolite <- metabolite
        df <- df[df$P.adj < 0.05, ]
        if (nrow(df) > 0) return(df)
        else return(NULL)
    }))
    
    # 4. Merge with metabolite descriptions (if available)
    if (!is.null(dunn_sig)) {
        dunn_sig <- data.frame(dunn_sig, row.names = NULL)
        
        # Ensure metabo_descp is a tibble with a "Metabolite" column
        if (!("Metabolite" %in% colnames(metabo_descp))) {
            stop("metabo_descp must be a tibble/dataframe with a 'Metabolite' column.")
        }
        
        dunn_sig_descp <- dunn_sig %>%
            left_join(metabo_descp, by = "Metabolite")
        
        return(dunn_sig_descp)
        
    } else {
        message("No significant pairwise differences found after adjustment.")
        return(NULL)
    }
}

# ---------------------------------------------------------------
# Example Usage: Direct comparison - dunnTest
# ---------------------------------------------------------------
pos_dunn <- dunnTest_metabo(pos_krus, log2_matrix_pos, Metadata_55, pos_descp, "pos")
neg_dunn <- dunnTest_metabo(neg_krus, log2_matrix_neg, Metadata_55, neg_descp, "neg")
