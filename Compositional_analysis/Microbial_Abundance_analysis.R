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


# ------------------------------
# 1. Set Working Directory
# ------------------------------

setwd("~/Documents/Microbiome_Analysis_R/4.Compositional_Profiles/Stat")

# ---- Load and Prepare Data ----
taxa <- asv_genus_with_metadata_rel %>% select(-SampleID, -Group)  # example; repeat for each level
sample_metadata$Group <- factor(sample_metadata$Group, levels = c("N", "OW", "OB"))

# ---- Kruskal-Wallis Test ----
krus <- lapply(taxa, function(x) kruskal.test(x ~ Group, data = sample_metadata))
krus_pvalue <- data.frame(
  Taxa = colnames(taxa),
  p.value = sapply(krus, getElement, name = "p.value")
) %>% mutate_if(is.numeric, round, digit = 6)
Export(krus_pvalue, "5_krus_genus_pval.txt")

# ---- Dunn’s Post Hoc Test ----
library(FSA)
krus_sig_groups <- krus_pvalue %>% filter(p.value < 0.05)
dunn_test <- setNames(
  lapply(krus_sig_groups$Taxa, function(group) {
    dunn_res <- dunnTest(taxa[[group]] ~ Group, data = sample_metadata, method = "bh")
    as.data.frame(dunn_res$res)
  }),
  krus_sig_groups$Taxa
)
dunn_sig <- do.call(rbind, lapply(names(dunn_test), function(Taxa) {
  df <- dunn_test[[Taxa]]
  df$Taxa <- Taxa
  df[df$P.adj < 0.05, ]
}))
dunn_sig_bac <- data.frame(dunn_sig, row.names = NULL)
Export(dunn_sig_bac, "5_dunn_sig_genus.txt")

# ---- Boxplot Visualization ----
library(reshape2)
library(ggpubr)
abundance_long <- melt(asv_genus_with_metadata_rel, id.vars = c("SampleID", "Group"),
                       variable.name = "Taxa", value.name = "Abundance")
sig_taxa <- unique(dunn_sig_bac$Taxa)
abundance_sig <- abundance_long %>%
  filter(Taxa %in% sig_taxa, Taxa != "unclassified") %>%
  mutate(Group = factor(Group, levels = c("N", "OW", "OB")))

ggplot(abundance_sig, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.7, alpha = 0.6) +
  facet_wrap(~ Taxa, scales = "free_y") +
  scale_fill_manual(values = c("N" = "grey", "OW" = "#FFA500", "OB" = "darkred")) +
  labs(y = "Relative Abundance", x = "Group") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold.italic", size = 9),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.grid = element_blank())

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
