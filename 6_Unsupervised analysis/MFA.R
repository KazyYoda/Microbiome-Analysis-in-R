################################################################################
# Microbiome Analysis: Gut Microbiota + Metabolites - Multiple Factor Analysis
# Author: Lucky
# Date: 2025-07
# Description:
#   1. Prepare CLR-transformed genus data and scaled metabolite data (pos/neg ions)
#   2. Merge datasets by sample ID and format for MFA input
#   3. Perform MFA using the FactoMineR package
#   4. Visualize sample projections on MFA dimensions (Dim 1 & 2)
#   5. Extract and plot variable contributions by data group (microbiota vs. metabolites)
#   6. Identify top contributing features and export their coordinates
################################################################################


# Set working directory
setwd("~/Documents/Obese_Microbiome/7.MFA")

# Load sPLS-DA workspace
load("~/Documents/Obese_Microbiome/6.sPLS-DA/sPLSDA.RData")


#---------------------------------------------------------
# 1. 🧹 Data Preparation for MFA
#---------------------------------------------------------

# Load required libraries
library(dplyr)
library(readxl)
library(compositions)
library(FactoMineR)
library(factoextra)
library(ggplot2)

# 📑 Load genus count and Kruskal-Wallis results
Counts_Genus <- read_excel("5_Counts_Genus.xlsx")
krus_genus_pval <- read.delim("5_krus_genus_pval.txt")

# 🎯 Filter significant genera (p < 0.05)
sig_kruskal_genus <- krus_genus_pval %>%
  filter(p.value < 0.05) %>%
  pull(Taxa)

# 🧬 Subset genus count data
sig_genus <- Counts_Genus %>%
  select(1:2, all_of(sig_kruskal_genus))

# 🔗 Match SampleID and Group to metadata
matched_sig_genus_55 <- Metadata_55 %>%
  left_join(sig_genus, by = c("SampleID", "Group"))

# ✅ Confirm SampleID alignment
stopifnot(identical(matched_sig_genus_55$SampleID, Metadata_55$SampleID))

# 🧮 Remove metadata columns & apply CLR transformation
genus_counts_pseudo <- matched_sig_genus_55 %>%
  select(-SampleID, -Group) + 1

genus_clr <- clr(genus_counts_pseudo)

# 🧪 Combine metabolomics + genus CLR data
metabo_full <- cbind(data.frame(Group = Metadata_55$Group,
                                metabo_scaled, 
                                genus_clr))

# ⚙️ Reorder factor levels
metabo_full$Group <- factor(metabo_full$Group, levels = c("N", "OW", "OB"))


#---------------------------------------------------------
# 2. 📊 Perform MFA
#---------------------------------------------------------

# MFA input structure
res.mfa <- MFA(
  metabo_full,
  group = c(1, 95, 38, 26),
  type = c("n", "c", "c", "c"),
  name.group = c("Group", "pos_ion", "neg_ion", "Genus"),
  num.group.sup = 1,
  graph = FALSE
)

# 📈 Scree plot
fviz_screeplot(res.mfa)

#---------------------------------------------------------
# 3. 🔍 Group Variable Contributions
#---------------------------------------------------------

# Group coordinates and contributions
group <- get_mfa_var(res.mfa, "group")
head(group$coord)
head(group$contrib)

# 🧭 Variable group map (Dim 1 & 2)
fviz_mfa_var(res.mfa, "group", axes = c(1,2), repel = TRUE) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")

# 🧭 Variable group map (Dim 3 & 4)
fviz_mfa_var(res.mfa, "group", axes = c(3,4), repel = TRUE) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")

# 📊 Group contribution barplots
fviz_contrib(res.mfa, "group", axes = 1) + theme_bw()
fviz_contrib(res.mfa, "group", axes = 2) + theme_bw()
fviz_contrib(res.mfa, "group", axes = 3) + theme_bw()


#---------------------------------------------------------
# 4. 👥 Individual Sample Projection
#---------------------------------------------------------

group_colors <- c(N = "grey", OW = "#FFA500", OB = "darkred")

fviz_mfa_ind(res.mfa, axes = c(1,2),
             label = "none",
             habillage = "Group",
             palette = group_colors,
             addEllipses = TRUE,
             ellipse.type = "confidence",
             repel = TRUE) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12))


#---------------------------------------------------------
# 5. 🧪 Quantitative Variable Analysis
#---------------------------------------------------------

# 📌 Variable correlation map (Dim 1 & 2)
fviz_mfa_var(res.mfa, "quanti.var", axes = c(1,2), 
             palette = "jco", repel = TRUE, geom = "point") +
  theme_bw()

# 📌 Variable correlation map (Dim 3 & 4)
fviz_mfa_var(res.mfa, "quanti.var", axes = c(3,4), 
             palette = "jco", repel = TRUE, geom = "point") +
  theme_bw()

# 📊 Contributions to dimensions (Top 50 variables)
for (dim in 1:4) {
  fviz_contrib(res.mfa, choice = "quanti.var", axes = dim, top = 50,
               palette = "jco") +
    theme_bw() +
    theme(axis.text = element_text(size = 8, angle = 45, hjust = 1),
          legend.text = element_text(size = 9))
}

#---------------------------------------------------------
# 6. 📌 Dimension Descriptions
#---------------------------------------------------------
dimdesc(res.mfa, axes = 1)
dimdesc(res.mfa, axes = 2)
dimdesc(res.mfa, axes = 3)
dimdesc(res.mfa, axes = 4)
