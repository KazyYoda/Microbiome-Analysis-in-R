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

setwd("~/Documents/Microbiome_Analysis_R/4.Compositional_Profiles")
load("~/Documents/Microbiome_Analysis_R/4.Compositional_Profiles/Compositional_Profiles.RData")


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

  # ---- Step 4: Boxplot for Significant Taxa ----
  abundance_long <- melt(taxa_raw, id.vars = c("SampleID", "Group"),
                         variable.name = "Taxa", value.name = "Abundance")

  sig_taxa <- unique(dunn_sig$Taxa)

  abundance_sig <- abundance_long %>%
    filter(Taxa %in% sig_taxa, Taxa != "unclassified") %>%
    mutate(Group = factor(Group, levels = c("N", "OW", "OB")))

  # Plot boxplots for significant taxa
  p <- ggplot(abundance_sig, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.7, alpha = 0.6) +
    facet_wrap(~ Taxa, scales = "free_y") +
    scale_fill_manual(values = c("N" = "grey", "OW" = "#FFA500", "OB" = "darkred")) +
    labs(title = paste("Differential Abundance -", level_name),
         y = "Relative Abundance", x = "Group") +
    theme_bw() +
    theme(strip.text = element_text(face = "bold.italic", size = 9),
          legend.position = "none",
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 9),
          panel.grid = element_blank())
    
    print(p)
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

               
# ---- Heatmap Visualization ----
library(ComplexHeatmap)
library(circlize)
group_colors <- c(N = "grey", OW = "#FFA500", OB = "darkred")
Rel_Genus_sig <- asv_genus_with_metadata_rel %>% select(SampleID, Group, all_of(sig_taxa))
Rel_Genus_sig_data <- Rel_Genus_sig %>% select(-SampleID, -Group) %>% as.matrix()
Rel_Genus_sig_data_t <- t(Rel_Genus_sig_data)
colnames(Rel_Genus_sig_data_t) <- Rel_Genus_sig$SampleID
sample_group <- data.frame(Group = sample_metadata$Group)
rownames(sample_group) <- sample_metadata$SampleID
ha <- HeatmapAnnotation(df = sample_group, col = list(Group = group_colors))
col_fun <- colorRamp2(c(min(Rel_Genus_sig_data_t), mean(Rel_Genus_sig_data_t), max(Rel_Genus_sig_data_t)),
                      c("#89ABE3", "#FDFD96", "#FC766A"))
Heatmap(Rel_Genus_sig_data_t, name = "Relative abundance", top_annotation = ha,
        col = col_fun, border = TRUE, cluster_rows = TRUE, cluster_columns = TRUE,
        column_split = sample_group$Group,
        row_names_gp = gpar(fontsize = 8, fontface = "italic"),
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(direction = "vertical"))
