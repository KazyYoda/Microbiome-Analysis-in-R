################################################################################
# Microbiome Analysis: Gut Microbiota + Metabolites - Redundancy Analysis (RDA)
# Author: Lucky
# Date: 2025-07
# Description:
#   1. Prepare genus-level microbiota and scaled metabolite datasets
#   2. Apply CLR transformation to genus data (after pseudocount)
#   3. Perform RDA with BMI group as the explanatory variable
#   4. Test model and axis significance via permutation-based ANOVA
#   5. Visualize top contributing features using a clean ggplot2 biplot (Scaling 2)
#   6. Highlight genera and metabolites (positive/negative ions) with arrows
################################################################################

#-------------------------------
# 0. Setup Environment
#-------------------------------

# Set working directory
setwd("~/Documents/Microbiome_Analysis_R/6.RDA")

# Load processed metabolite and metadata objects
load("~/Documents/Microbiome_Analysis_R/5.sPLS-DA/sPLSDA.RData")

# Load required packages
library(dplyr)
library(tidyr)
library(readxl)
library(compositions)
library(vegan)
library(ggplot2)
library(ggrepel)
library(permute)
library(lattice)

#-------------------------------
# 1. Prepare Genus Data (CLR)
#-------------------------------

# Load genus counts
Counts_Genus <- read_excel("Documents/Microbiome_Analysis_R/3.Compositional_Profiles/5_Counts_Genus.xlsx")

# Filter to 55 samples used in metabolite analysis
Genus <- Counts_Genus %>% filter(SampleID %in% Metadata_55$SampleID)

# Join metadata to ensure consistent sample order
matched_genus_55 <- Metadata_55 %>%
  left_join(Genus, by = c("SampleID", "Group"))

# Confirm alignment of sample order
identical(matched_genus_55$SampleID, Metadata_55$SampleID)

# Remove metadata columns and add pseudocount
Genus_55 <- matched_genus_55 %>%
  select(-SampleID, -Group) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0) + 1))

# CLR transformation (requires numeric matrix)
# ✅ genus = apply CLR transformation
genus_clr <- clr(Genus_55)




#-------------------------------
# 2. Redundancy Analysis (RDA)
#-------------------------------

#📍 ----- RDA on metabolites -----
# •	✅ Response matrix: metabo_scaled
# •	✅ Explanatory variable: Metadata$Group (BMI Group)

# Define RDA model
# ✅ metabo_scaled = log-transformed + scaled
rda_model <- rda(metabo_scaled ~ Group, data = Metadata_55)


#📍 ----- RDA with multi-omics as response -----
# Associations between genus and metabolites across BMI groups
# •	✅ Response matrix: both genus_clr and metabo_scaled
# •	✅ Explanatory variable: Metadata$Group (e.g. BMI group)

# Combine CLR-genus and scaled metabolites
response_matrix <- cbind(genus_clr, metabo_scaled)

# Define RDA model
# ✅ genus_clr = clr-transformed 
# ✅ metabo_scaled = log-transformed + scaled
rda_model <- rda(response_matrix ~ Group, data = Metadata_55)


# Model summary and adjusted R²
summary(rda_model)
# • Constrained Proportion: variance of Y explained by X 
# • Unconstrained Proportion: unexplained variance in Y 

# Canonical coefficients
coef(rda_model)
RsquareAdj(rda_model)


# Permutation tests for significance (PERMUTATION-BASED ANOVA)
# Test overall model
set.seed(123)
anova.cca(rda_model, permutations = 999)

# Test each constrained axis
set.seed(123)
anova.cca(rda_model, by = "axis", permutations = 999)

# Test terms (Group, metabolites, etc.)
set.seed(123)
anova.cca(rda_model, by = "terms", permutations = 999)

# Quick base RDA plot
plot(rda_model, scaling = 1)

#🔬 Interpretation

# • Scaling = 1 = distances between samples matter
# •	Scaling = 2 = angles between variables (species) matter (better when you want to explore relationships between features like genus or metabolites).

#This RDA will help you answer:
#•	✅ How much variance in combined genus + metabolite profiles is explained by BMI group?
#•	✅ Are there visible clusters or shifts in profiles across BMI categories?
#•	✅ Which genera and metabolites drive the separation?




#-------------------------------
# 3. Visualization - ggplot2
#-------------------------------

# ------- Helper function ------
rda_plot <- function(rda_model, metadata, top_n = 100, scale_factor = 2, scaling_n = 2) {
  
  # Extract RDA scores with desired scaling
  rda_scores <- scores(rda_model, display = c("sites", "species"), scaling = scaling_n)
  
  # Prepare site data (sample scores + group info)
  sites_df <- as.data.frame(rda_scores$sites) %>%
    mutate(Group = metadata$Group)
  
  # Extract species scores and compute arrow lengths
  species_scores <- as.data.frame(rda_scores$species)
  arrow_lengths <- sqrt(rowSums(species_scores^2))
  
  # Select top N species based on arrow length
  top_species <- names(sort(arrow_lengths, decreasing = TRUE))[1:top_n]
  species_df <- species_scores[top_species, ]
  
  # Scale species scores to adjust arrow length
  species_df_scaled <- species_df * scale_factor
  
  # Add feature classification: Genus, positive_ion, negative_ion
  species_df_scaled$ion_type <- case_when(
    grepl("^pos", rownames(species_df_scaled)) ~ "positive_ion",
    grepl("^neg", rownames(species_df_scaled)) ~ "negative_ion",
    TRUE ~ "Genus" # Default to "Genus" if it's neither pos nor neg
  )
  
  # Compute explained variance for axis labels
  rda_var <- summary(rda_model)$concont$importance["Proportion Explained", 1:2]
  
  # Define color palettes
  group_colors <- c(N = "grey", OW = "#FFA500", OB = "darkred")
  
  if (scaling_n == 1) {
    arrow_colors <- c("Genus" = "grey", "negative_ion" = "grey", "positive_ion" = "grey")
  } else {
    arrow_colors <- c("Genus" = "#660099", "negative_ion" = "#e02b35", "positive_ion" = "#228B22")
  }
  
  # Create the biplot
  p <- ggplot() +
    geom_point(data = sites_df, aes(RDA1, RDA2, color = Group), size = 2, alpha = 0.6) +
    geom_segment(data = species_df_scaled, 
                 aes(x = 0, y = 0, xend = RDA1, yend = RDA2, color = ion_type),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 0.5) +
    geom_text_repel(data = species_df_scaled, 
                    aes(RDA1, RDA2, label = rownames(species_df_scaled), color = ion_type),
                    size = 3, segment.colour = "grey", max.overlaps = 100) +
    scale_color_manual(values = c(arrow_colors, group_colors)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
    labs(
      title = paste("RDA Biplot (Scaling", scaling_n, "): Top", top_n, "Features"),
      x = paste0("RDA1 (", round(100 * rda_var[1], 2), "%)"),
      y = paste0("RDA2 (", round(100 * rda_var[2], 2), "%)")
    ) +
    theme_bw() +
    theme(legend.position = "right")
  
  return(p)
}


# ---------------------------------------------------------------
# Example Usage: RDA plot
# ---------------------------------------------------------------

# Parameters:
# • rda_model - RDA model object (from vegan::rda())
# • metadata - Metadata with SampleID and Group columns
# • top_n - Number of top species to display as arrows
# • scale_factor - Multiplier to stretch arrow length
# • scaling_n - 1 = distances matter (samples), 2 = angles matter (features)
# 📝 When scaling_n = 1, all feature arrows are grey to reflect distance-based scaling.

rda_plot(rda_model = rda_model, metadata = Metadata_55, top_n = 100, 
         scale_factor = 2, scaling_n = 2)

rda_plot(rda_model = rda_model, metadata = Metadata_55, top_n = 100, 
         scale_factor = 2, scaling_n = 1)

# Multiple plots
library(gridExtra)
grid.arrange(
  sc1, 
  sc2, 
  layout_matrix = cbind(c(1, 2))
)
