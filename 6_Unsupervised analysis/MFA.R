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
setwd("~/Documents/Microbiome_Analysis_R/7.MFA")


# Load sPLS-DA workspace
load("~/Documents/Microbiome_Analysis_R/5.sPLS-DA/sPLSDA.RData")


# Load required libraries
library(dplyr)
library(readxl)
library(compositions)
library(FactoMineR)
library(factoextra)
library(ggplot2)


#---------------------------------------------------------
# 1. Data Preparation for MFA
#---------------------------------------------------------

# Load genus count and Kruskal-Wallis results
Counts_Genus <- read_excel("Documents/Microbiome_Analysis_R/3.Compositional_Profiles/5_Counts_Genus.xlsx")
krus_genus_pval <- read.delim("~/Documents/Obese_Microbiome/3.Compositional_Profiles/Stat/5_krus_genus_pval.txt")


# Filter significant genera (p < 0.05)
sig_kruskal_genus <- krus_genus_pval %>%
  filter(p.value < 0.05) %>%
  pull(Taxa)


# Subset genus count data
sig_genus <- Counts_Genus %>%
  select(1:2, all_of(sig_kruskal_genus))


# Match SampleID and Group to metadata
matched_sig_genus_55 <- Metadata_55 %>%
  left_join(sig_genus, by = c("SampleID", "Group"))


# Confirm SampleID alignment
identical(matched_sig_genus_55$SampleID, Metadata_55$SampleID)


# Remove metadata columns & apply CLR transformation
genus_counts_pseudo <- matched_sig_genus_55 %>%
  select(-SampleID, -Group) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0) + 1))

genus_clr <- clr(genus_counts_pseudo)


# Combine metabolomics + genus CLR data
metabo_full <- cbind(data.frame(Group = Metadata_55$Group,
                                metabo_scaled, 
                                genus_clr))

# Reorder factor levels
metabo_full$Group <- factor(metabo_full$Group, levels = c("N", "OW", "OB"))




#---------------------------------------------------------
# 2. Perform MFA
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

# Scree plot
fviz_screeplot(res.mfa)




#---------------------------------------------------------
# 3. Group Variable Contributions
#---------------------------------------------------------

# Eigenvalues / Variances:
eig.val <- get_eigenvalue(res.mfa)
head(eig.val)


# Group coordinates and contributions
group <- get_mfa_var(res.mfa, "group")

# Coordinates of groups
head(group$coord)

# Cos2: quality of representation on the factor map
head(group$cos2)

# Contributions to the dimensions
head(group$contrib)


# Variable group map (Dim 1 & 2)
fviz_mfa_var(res.mfa, "group", axes = c(1,2), repel = TRUE) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")

# Variable group map (Dim 3 & 4)
fviz_mfa_var(res.mfa, "group", axes = c(3,4), repel = TRUE) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")


# Group contribution barplots
# To view Dim 2 and 3 manually, change "axes"
fviz_contrib(res.mfa, "group", axes = 1) + theme_bw() +
  theme(legend.position = "right", legend.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_blank())




#---------------------------------------------------------
# 4. Individual Sample Projection
#---------------------------------------------------------
ind <- get_mfa_ind(res.mfa)
ind

group_colors <- c(N = "grey", OW = "#FFA500", OB = "darkred")

# To view Dim 3 and 4 manually, change "axes = c(3,4)"
fviz_mfa_ind(res.mfa, axes = c(1,2), 
             labelsize = 2, 
             invisible = "none", #omit  the supplementary qualitative variable
             habillage = c("Group"), # colored by groups
             palette = group_colors, 
             label = "none",
             addEllipses = TRUE, 
             ellipse.type = "confidence", 
             repel = TRUE)+ # Avoid text overlapping)
  theme(aspect.ratio = 1,
        axis.title.y = element_text(face = "plain", size = rel(1.2),
                                    margin = margin(t = 0, r = 15, b = 0, l = 0)
                                    , hjust = 0.5), #label of y
        axis.title.x = element_text(face = "plain", size = rel(1.2),
                                    margin = margin(t = 15, r = 0, b = 0, l = 0)
                                    , vjust = 0.5),
        axis.text.y = element_text(face = "plain", size = rel(1.2), colour = "black"), #scale y label
        axis.text.x = element_text(face = "plain", size = rel(1.2), colour = "black",
                                   angle = 0,
                                   margin = margin(t = 0, r = 0, b = 0, l = 0),
                                   vjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.ticks.length.x = unit(0.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.spacing.y = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "cm"),
        legend.position = "right",
        legend.text = element_text(face = "plain", size = rel(1.2))
  )




#---------------------------------------------------------
# 5. Quantitative Variable Analysis
#---------------------------------------------------------

# Variable correlation map (Dim 1 & 2)
fviz_mfa_var(res.mfa, "quanti.var", axes = c(1,2), 
             palette = "jco", repel = TRUE, geom = "point") +
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14))


# Variable correlation map (Dim 3 & 4)
fviz_mfa_var(res.mfa, "quanti.var", axes = c(3,4), 
             palette = "jco", repel = TRUE, geom = "point") +
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14))


# Contributions to dimensions (Top 50 variables)
for (dim in 1:4) {
  p <- fviz_contrib(res.mfa, choice = "quanti.var", axes = dim, top = 50,
               palette = "jco") +
    theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 9, colour = "black"),
        axis.text = element_text(color = "black", size = 8, angle = 45, hjust = 1),
        axis.title.x = element_blank())
  print(p)
}


#---------------------------------------------------------
# 6. Dimension Descriptions
#---------------------------------------------------------
dimdesc(res.mfa, axes = 1)
dimdesc(res.mfa, axes = 2)
dimdesc(res.mfa, axes = 3)
dimdesc(res.mfa, axes = 4)
