############################################################
# Microbiome Analysis in R: Beta Diversity Analysis
# Author: Lucky
# Date: [2025-7]
# Description:
#   This script calculates beta diversity metrics, performs 
#   PCoA ordination, visualizes results, and runs PERMANOVA 
#   and beta-dispersion tests for group differences.
############################################################


# ------------------------------
# 1. Setup and Load Data
# ------------------------------
setwd("~/Documents/Microbiome_Analysis_R/3.Beta_diversity")
load("~/Documents/Microbiome_Analysis_R/1.Raw_Data/Building_Phyloseq.RData")


# Load packages
library(phyloseq)
library(vegan)
library(picante)
library(ggplot2)
library(pairwiseAdonis)
library(gridExtra)




# ------------------------------
# 2. Prepare Phyloseq Components
# ------------------------------
otu_table <- otu_table(ps)
tax_table <- tax_table(ps)
sample_data <- sample_data(ps)


sample_data$Group <- factor(sample_data$Group, levels = c("N", "OW", "OB"))
tree <- phy_tree(ps)


# Optional tree quality check
is.rooted(tree)
is.binary(tree)  # UniFrac assumes binary trees ((meaning some nodes have more than 2 descendants)




# ------------------------------
# 3. Compute Beta Diversity Distances
# ------------------------------
unw_unifrac <- distance(ps, method = "unifrac", weighted = FALSE)
w_unifrac <- distance(ps, method = "unifrac", weighted = TRUE)
bray_curtis <- distance(otu_table, method = "bray")




# ------------------------------
# 4. PCoA Ordination
# ------------------------------
unw_pcoa   <- ordinate(ps, method = "PCoA", distance = unw_unifrac)
w_pcoa     <- ordinate(ps, method = "PCoA", distance = w_unifrac)
bray_pcoa  <- ordinate(ps, method = "PCoA", distance = bray_curtis)




# ------------------------------
# 5. Visualization Function
# ------------------------------

# Function
plot_pcoa <- function(ord, ps, method_name, dist_name) {
  plot_ordination(ps, ord, color = "Group") +
    geom_point(size = 3) +
    stat_ellipse(aes(color = Group), level = 0.95) +
    scale_color_manual(values = c("N" = "grey", "OW" = "#FFA500", "OB" = "darkred")) +
    labs(
      title = paste("PCoA -", dist_name),
      x = paste0("PCoA1 (", round(100 * ord$values$Relative_eig[1], 2), "%)"),
      y = paste0("PCoA2 (", round(100 * ord$values$Relative_eig[2], 2), "%)")
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 10),
          legend.position = "none")
}

# Unweighted UniFrac plot
p1 <- plot_pcoa(unw_pcoa, ps, "PCoA", "Unweighted UniFrac")
p1

# Weighted UniFrac plot
p2 <- plot_pcoa(w_pcoa, ps, "PCoA", "Weighted UniFrac")
p2

# Bray-Curtis plot
p3 <- plot_pcoa(bray_pcoa, ps, "PCoA", "Bray-Curtis")
p3




# ---------------------------------------
# 6. Run PERMANOVA + pairwise PERMANOVA
# ---------------------------------------
# Extract the grouping variable directly
grouping <- factor(sample_data(ps)$Group,
                   levels = c("N", "OW", "OB"))


# -------------------------------------
# Run PERMANOVA and Pairwise PERMANOVA
# -------------------------------------

# Function
run_permanova <- function(dist_matrix, name, group_vector) {
  cat(paste0("\n=======================================\n"))
  cat(paste0(" PERMANOVA Results: ", name, "\n"))
  cat(paste0("=======================================\n"))

  # Set seed for reproducibility
  set.seed(123)

  # Main PERMANOVA test
  permanova <- adonis2(dist_matrix ~ group_vector, permutations = 999)
  print(permanova)

  # Pairwise PERMANOVA
  cat("\nPairwise PERMANOVA:\n")
  pairwise <- pairwise.adonis(dist_matrix, 
                               factors = group_vector, 
                               perm = 999, 
                               p.adjust.m = "BH")
  print(pairwise)
}


# Unweighted UniFrac 
run_permanova(unw_unifrac, "Unweighted UniFrac", grouping)

# Weighted UniFrac 
run_permanova(w_unifrac, "Weighted UniFrac", grouping)

# Bray-Curtis plot
run_permanova(bray_curtis, "Bray-Curtis", grouping)




# -----------------------------------
# 7. Beta Dispersion & Visualization
# -----------------------------------

# Function
plot_betadisper <- function(dist, grouping, method, sqrt = FALSE) {
  
  # Perform Beta-Dispersion (betadisper)
  disp <- betadisper(dist, grouping, sqrt.dist = sqrt)
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Perform statistical test for dispersion
  print(permutest(disp, permutations = 999))
  
  # Visualization
  df <- data.frame(Group = disp$group, 
                   Distance = disp$distances)
  
  p <- ggplot(df, aes(x = Group, y = Distance, fill = Group)) +
    geom_boxplot() +
    labs(title = paste("Beta Dispersion -", method),
         y = "Distance to centroid") +
    scale_fill_manual(name = "Group",
                       values = c("0MC" = "grey", "0MT" = "darkgoldenrod1", 
                                  "3MC" = "gray20", "3MT" = "darkorange2"))+
    theme_bw() +
    theme(panel.grid = element_blank(),
      legend.position = "right")
  
  return(p)
}


# Unweighted UniFrac 
d1 <- plot_betadisper(unw_unifrac, "Unweighted UniFrac")

# Weighted UniFrac
d2 <- plot_betadisper(w_unifrac, "Weighted UniFrac", sqrt = TRUE)

# sqrt transformation for UniFrac distance matrix that is non-Euclidean (some eigenvalues become negative)
# sqrt.dist = TRUE: applies a square root transformation to the distance matrix inside the betadisper() function
# This is useful for non-Euclidean distances, such as weighted UniFrac, to reduce negative eigenvalues.


# Bray-Curtis plot
d3 <- plot_betadisper(bray_curtis, "Bray-Curtis")




# ------------------------------
# 8. Combine All Plots
# ------------------------------
grid.arrange(p1, d1, 
             p2, d2, 
             p3, d3, 
             ncol = 2, nrow = 3)
