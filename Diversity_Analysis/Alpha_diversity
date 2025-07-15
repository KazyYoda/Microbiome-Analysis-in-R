############################################################
# Microbiome Analysis in R: Alpha Diversity Analysis
# Author: Lucky
# Date: [2025-7]
# Description:
#   This script computes alpha diversity metrics from a 
#   phyloseq object, performs statistical tests (Kruskal-Wallis), 
#   and generates boxplots for visualization.
############################################################


# ------------------------------
# 1. Set Working Directory
# ------------------------------
setwd("~/Documents/Obese_Microbiome/2.Alpha_diversity")


# ------------------------------
# 2. Load Phyloseq Object & Packages
# ------------------------------
load("~/Documents/Obese_Microbiome/1.Raw_Data/Building_Phyloseq.RData")

library(phyloseq)
library(car)
library(dplyr)


# ------------------------------
# 3. Estimate Alpha Diversity
# ------------------------------

# Standard metrics: Observed, Shannon, Chao1
alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon", "Chao1"))
head(alpha_div)

# Extract metadata
meta <- as(sample_data(ps), "data.frame")
meta$Group <- factor(meta$Group, levels = c("N", "OW", "OB"))

# Combine alpha diversity with metadata
alpha_div_meta <- cbind(meta, alpha_div)
head(alpha_div_meta)


# ------------------------------
# 4. Calculate Faith’s Phylogenetic Diversity (PD)
# ------------------------------
install.packages("picante")  # If not already installed
library(picante)

# Extract OTU table and tree from phyloseq
otu <- as(otu_table(ps), "matrix")
tree <- phy_tree(ps)
tree$tip.label <- as.character(tree$tip.label)

# Align OTU table rows with tree tip labels
otu <- otu[tree$tip.label, ]
otu_t <- t(otu)  # Transpose so rows = samples, columns = ASVs

# Compute PD
pd_result <- pd(otu_t, tree)
pd_result$SampleID <- rownames(pd_result)

# Merge PD with alpha diversity table
alpha_div_meta_combined <- merge(alpha_div_meta,
                                 pd_result[, c("SampleID", "PD")],
                                 by = "SampleID")

head(alpha_div_meta_combined)

# Export to file
Export(alpha_div_meta_combined, "Alpha_diversity.txt")


# ------------------------------
# 5. Statistical Testing (Kruskal-Wallis)
# ------------------------------

# Clean input: remove standard error of Chao1
alpha_df <- alpha_div_meta_combined %>% select(-se.chao1)
alpha_metrics <- c("Observed", "Chao1", "Shannon", "PD")

# Apply Kruskal-Wallis test
kruskal_results <- lapply(alpha_metrics, function(metric) {
  kruskal.test(as.formula(paste(metric, "~ Group")), data = alpha_df)
})
names(kruskal_results) <- alpha_metrics

# Extract p-values
kruskal_summary <- data.frame(
  Metric = alpha_metrics,
  P_Value = sapply(kruskal_results, function(x) x$p.value)
)

print(kruskal_summary)
Export(kruskal_summary, "kruskal_pval.txt")


# ------------------------------
# 6. Visualization: Boxplots
# ------------------------------
library(tidyr)
library(ggplot2)

# Convert to long format
alpha_long <- pivot_longer(alpha_df, 
                           cols = all_of(alpha_metrics),
                           names_to = "Index",
                           values_to = "Value")

# Order diversity indices
alpha_long$Index <- factor(alpha_long$Index, 
                           levels = c("Observed", "Chao1", "Shannon", "PD"))

# Boxplot
ggplot(alpha_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  facet_wrap(~ Index, scales = "free_y", nrow = 1) +
  labs(x = "Group", y = "Alpha Diversity Index") +
  scale_fill_manual(values = c("N" = "grey", "OW" = "#FFA500", "OB" = "darkred")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
