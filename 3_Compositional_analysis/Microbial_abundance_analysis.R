###########################################################################
# Microbiome Analysis in R: Microbial Abundance analysis
# Author: Lucky
# Date: [2025-07]
# Description:
#   1. Prepare relative abundance data by taxonomic rank
#   2. Perform Kruskal-Wallis test
#   3. Conduct Dunn’s post hoc test for significant taxa
#   4. Visualize results using boxplots, heatmaps, and stacked bar plots
##########################################################################


# ----------------------------------------
# 1. Set Working Directory & Load .RData
# ----------------------------------------

setwd("~/Documents/Microbiome_Analysis_R/3.Compositional_Profiles")
load("~/Documents/Microbiome_Analysis_R/3.Compositional_Profiles/Compositional_Profiles.RData")


# Load packages
library(readxl)
library(dplyr)
library(FSA)
library(reshape2)
library(ggpubr)
library(rio)
library(car)




# ---------------------------------------------------------------
# 2. Abundance Analysis and Visualization by Taxonomic Rank
# ---------------------------------------------------------------

# Set working directory
setwd("~/Documents/Obese_Microbiome/4.Compositional_Profiles/Stat")


# Use sample metadata from "Compositional_Profiles.RData"
sample_metadata <- sample_metadata %>%
  mutate(Group = factor(Group, levels = c("N", "OW", "OB")))


# Helper function for statistical analysis and visualization
stat_taxonomic_level <- function(level_name, prefix, input_file) {
  message("▶️ Processing: ", level_name)

  # ---- Step 1: Load and Prepare Data ----
  taxa_raw <- read_excel(input_file)
  taxa <- taxa_raw %>% select(-SampleID, -Group)

  # ---- Step 2: Kruskal-Wallis Test ----
  krus <- lapply(taxa, function(x) kruskal.test(x ~ sample_metadata$Group))
  krus_pvalue <- data.frame(
    Taxa = colnames(taxa),
    p.value = sapply(krus, function(x) x$p.value)
  ) %>%
    mutate(p.value = round(p.value, 6))

  # Export p-values 
  Export(krus_pvalue, paste0(prefix, "_krus_", level_name, "_pval.txt"))

  # ---- Step 3: Dunn’s Post Hoc Test (BH adjusted) ----
  krus_sig <- krus_pvalue %>% filter(p.value < 0.05)

  dunn_test <- setNames(
    lapply(krus_sig$Taxa, function(taxa_name) {
      dunn_res <- dunnTest(taxa[[taxa_name]] ~ sample_metadata$Group, method = "bh")
      df <- as.data.frame(dunn_res$res) # Convert results to dataframe
      df$Taxa <- taxa_name
      df
    }),
    krus_sig$Taxa # Assign names based on taxa
  )

   # Filter significant pairs (adjusted p-value < 0.05)
  dunn_sig <- bind_rows(
    lapply(dunn_test, function(df) df %>% filter(P.adj < 0.05))
  )

  # Export Dunn results 
  Export(dunn_sig, paste0(prefix, "_dunn_sig_", level_name, ".txt"))

  # ---- Step 4: Conditional Boxplot ----
  if (level_name != "Genus") {
    abundance_long <- melt(taxa_raw, id.vars = c("SampleID", "Group"),
                           variable.name = "Taxa", value.name = "Abundance")

    sig_taxa <- unique(dunn_sig$Taxa)

    abundance_sig <- abundance_long %>%
      filter(Taxa %in% sig_taxa, Taxa != "unclassified") %>%
      mutate(Group = factor(Group, levels = c("N", "OW", "OB")))

    p <- ggplot(abundance_sig, aes(x = Group, y = Abundance, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 0.7, alpha = 0.6) +
      facet_wrap(~ Taxa, scales = "free_y") +
      scale_fill_manual(values = c("N" = "grey", "OW" = "#FFA500", "OB" = "darkred")) +
      labs(title = paste("Significant difference in abundance -", level_name),
           y = "Relative Abundance", x = "Group") +
      theme_bw() +
      theme(strip.text = element_text(face = "bold.italic", size = 8),
            legend.position = "none",
            axis.text = element_text(size = 9),
            axis.title = element_text(size = 9),
            panel.grid = element_blank())

    print(p)
  } else {
    message("⏭️ Skipping boxplot for Genus level.")
  }
}


# -----------------------------
# Run for all taxonomic levels
# -----------------------------
levels <- c("Phylum", "Class", "Order", "Family", "Genus")
prefixes <- c("1", "2", "3", "4", "5")
input_files <- c("1_Rel_Phylum.xlsx", "2_Rel_Class.xlsx", "3_Rel_Order.xlsx",
                 "4_Rel_Family.xlsx", "5_Rel_Genus.xlsx")

# Loop through levels
for (i in seq_along(levels)) {
  stat_taxonomic_level(level_name = levels[i],
                       prefix = prefixes[i],
                       input_file = input_files[i])
}

               
# ---- Heatmap Visualization of Significant Genus-Level Abundance ----
library(ComplexHeatmap)
library(circlize)


# ---- Step 1: Load Data ----
# Relative abundance table at Genus level (please check your directory)
Genus <- read_excel("5_Rel_Genus.xlsx")

    
# List of significantly different genera from Dunn's post hoc test
sig_taxa <- read.delim("5_dunn_sig_Genus.txt")

    
# ---- Step 2: Subset Data for Significant Genera ----
Rel_Genus_sig <- Genus %>%
  select(SampleID, Group, all_of(sig_taxa$Taxa))  # Keep only significant taxa

    
# Create a matrix for heatmap (Taxa as rows, Samples as columns)
Rel_Genus_sig_data <- Rel_Genus_sig %>%
  select(-SampleID, -Group) %>%
  as.matrix()

# Transpose the matrix (rows: taxa, columns: samples)
Rel_Genus_sig_data_t <- t(Rel_Genus_sig_data)
colnames(Rel_Genus_sig_data_t) <- Rel_Genus_sig$SampleID

    
# ---- Step 3: Prepare Group Annotation ----
group_colors <- c(N = "grey", OW = "#FFA500", OB = "darkred")

    
# Match group info by sample
sample_group <- data.frame(Group = sample_metadata$Group)
rownames(sample_group) <- sample_metadata$SampleID

    
# Annotation for sample groups
ha <- HeatmapAnnotation(df = sample_group, col = list(Group = group_colors))

    
# ---- Step 4: Define Heatmap Color Scale ----
col_fun <- colorRamp2(
  c(min(Rel_Genus_sig_data_t), 
    mean(Rel_Genus_sig_data_t), 
    max(Rel_Genus_sig_data_t)),
  c("#89ABE3", "#FDFD96", "#FC766A")
)

    
# ---- Step 5: Generate Heatmap ----
Heatmap(
  Rel_Genus_sig_data_t,
  name = "Relative abundance",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = col_fun,
  top_annotation = ha,
  column_split = sample_group$Group,  # Split samples by group
  border = TRUE,
  border_gp = gpar(col = "grey"),
  rect_gp = gpar(col = "grey30", lwd = 0.5),
  column_names_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  clustering_distance_rows = "spearman",
  clustering_distance_columns = "spearman",
  heatmap_legend_param = list(direction = "vertical")
)


    

# ---------------------------------------------------------------
# Stacked Barplot of Microbial Relative Abundance 
# ---------------------------------------------------------------

# Define custom color palette (at least as many as your top taxa)
custom_colors <- c(
  "#082a54",  # Dark Blue
  "#e02b35",  # Red
  "#f0c571",  # Gold
  "#59a89c",  # Teal
  "#a559aa",  # Purple
  "#3f4c6b",  # Dark Blue variant
  "#d92e2f",  # Deep Red
  "#f2b94f",  # Bright Gold
  "#5c8f8f",  # Teal variant
  "#9c4a9d",  # Purple variant
  "#4b8fbe",  # Steel Blue
  "#d66f51",  # Salmon Red
  "#ffb643",  # Lemon Yellow
  "#72b7b1",  # Pastel Teal
  "#b15598",  # Lavender Purple
  "#2c3e6b",  # Navy Blue
  "#e64d58",  # Coral Red
  "#d9c258",  # Soft Gold
  "#61b9a4",  # Mint Teal
  "#aa478c",  # Plum Purple
  "#1f2d49",  # Deep Navy
  "#fd7c4d",  # Tangerine Red
  "#fcd703",  # Vibrant Yellow
  "#72b5b7",  # Ocean Teal
  "#8a5398",  # Magenta Purple
  "#cecece"   # Gray (last color)
)

  
# Custom theme
theme_stack <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm")
  )

    
# ---- Helper Function for Plotting ----
plot_stacked_bar <- function(data, group_order, top_n = 10, level = "Phylum") {
  
  data$Group <- factor(data$Group, levels = group_order)

  # Step 1: Compute top N taxa by overall mean
  top_taxa <- colnames(data)[3:(2 + top_n)]

  # Step 2: Melt data to long format and collapse others
  data_long <- melt(data, id.vars = c("SampleID", "Group"),
                    variable.name = "Taxa", value.name = "Abundance") %>%
    mutate(Taxa = ifelse(as.character(Taxa) %in% top_taxa, as.character(Taxa), "Others")) 

  # Step 3: Compute mean abundance by group and taxa
  mean_abundance <- data_long %>%
    group_by(Group, Taxa) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop")

  # Taxa in reverse order 
  order_taxa <- c("Others", rev(top_taxa))
  
  # Step 4: Plot
  p <- ggplot(mean_abundance, aes(x = Group, y = MeanAbundance, 
                                  fill = factor(Taxa, level = order_taxa))) +
    geom_bar(stat = "identity", position = "fill") +
    #scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +  # manual percent formatting
    scale_fill_manual(values = custom_colors) +
    labs(x = "", y = "Relative Abundance", fill = level) +
    theme_stack

  return(p)
}

    
# ---------------------------------------------------------------
# Example Usage:
# ---------------------------------------------------------------

# --- Phylum level ---                   
Rel_Phylum <- read_excel("1_Rel_Phylum.xlsx") # please check your directory
phylum_plot <- plot_stacked_bar(data = Rel_Phylum, 
                                group_order = c("N", "OW", "OB"), 
                                top_n = 10, level = "Phylum")
print(phylum_plot)

                       
# --- Class level ---
Rel_Class <- read_excel("2_Rel_Class.xlsx") # please check your directory
class_plot <- plot_stacked_bar(data = Rel_Class, 
                                group_order = c("N", "OW", "OB"), 
                                top_n = 10, level = "Class")
print(class_plot)
