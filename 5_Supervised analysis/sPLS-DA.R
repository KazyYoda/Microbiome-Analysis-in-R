##################################################################################
# Microbiome Analysis in R: Metabolite Classification (Positive & Negative Ions)
# Author: Lucky
# Date: 2025-7
# Description:
#  This pipeline performs supervised multivariate analysis using sparse Partial 
#  Least Squares Discriminant Analysis (sPLS-DA) to identify metabolite signatures
#  that differentiate between sample groups.
###################################################################################


# Set Working Directory
setwd("~/Documents/Microbiome_Analysis_R/5.sPLS-DA")

# Load Required Libraries
library(readxl)
library(mixOmics)
library(MASS)
library(lattice)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Load Data
load("~/Documents/Microbiome_Analysis_R/5.Metabolites/Metabolites.RData")


#------------------------------------------------------
# 1. Data Preparation
#------------------------------------------------------
Merged_Metabo <- cbind(posmetabo, negmetabo)
Y <- Metadata_55$Group  # Must be a factor

# Log2 transform (handle zeroes beforehand if needed)
log_metabo <- log2(Merged_Metabo)
metabo_scaled <- scale(log_metabo)


#------------------------------------------------------
# 2. Initial sPLS-DA Model (2 components)
#------------------------------------------------------
splsda_result <- splsda(metabo_scaled, Y, ncomp = 2)

# Plot Individuals
plotIndiv(splsda_result, comp = 1:2, group = Y, legend = TRUE)

# Plot Variables
plotVar(splsda_result, comp = 1:2)

# VIP Scores
vip_scores <- vip(splsda_result)
print(vip_scores)


#------------------------------------------------------
# 3. Model Tuning
#------------------------------------------------------

# a) Tune Number of Components
perf.plsda <- perf(
  splsda_result, validation = 'Mfold',
  folds = 5, nrepeat = 100, auc = TRUE,
  progressBar = FALSE
)

# Extract the optimal component number
optimal.ncomp <- perf.plsda$choice.ncomp["BER", "max.dist"]
plot(perf.plsda, overlay = 'measure', sd = TRUE)


# b) Tune keepX
set.seed(123)
tuned <- tune.splsda(
  metabo_scaled, Y, ncomp = 2, validation = 'Mfold', folds = 5,
  dist = 'max.dist', measure = 'BER', nrepeat = 100,
  test.keepX = 1:133, progressBar = TRUE
)

# Extract the optimal component number and optimal feature count per component
optimal.keepX <- tuned$choice.keepX[1:2]
optimal.keepX

optimal.ncomp = tuned$choice.ncomp$ncomp 
optimal.ncomp

plot(tuned)




#------------------------------------------------------
# 4. Final Model & Visualization
#------------------------------------------------------
final_model <- splsda(metabo_scaled, Y, ncomp = 2, keepX = optimal.keepX)


group_colors <- c(N = "grey", OW = "#FFA500", OB = "darkred")

# Get plotIndiv output (coordinates and grouping)
plotIndiv(final_model,
          comp = c(1, 2),
          group = Y,
          col.per.group = group_colors,
          ind.names = FALSE,
          ellipse = FALSE, 
          legend = TRUE,
          legend.title = "Group",
          title = 'Fecal Metabolites, sPLS-DA Comps 1 & 2')

# ggplot2: Score Plot with Ellipses
indiv_df <- data.frame(
  comp1 = final_model$variates$X[, 1],
  comp2 = final_model$variates$X[, 2],
  Group = Y
)

# Obtain explained variance from "plotIndiv"
indv <- ggplot(indiv_df, aes(x = comp1, y = comp2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Group), linewidth = 1) +
  scale_color_manual(values = group_colors) +
  labs(title = "Fecal metabolites, sPLS-DA Comps 1 & 2",
       x = "Component 1 (X% expl. var)", y = "Component 2 (Y% explained. var)") +
  theme_bw(base_size = 8)


# Loading Plots
plotLoadings(final_model, comp = 1, method = 'mean', contrib = 'max', 
             size.name = 0.8, title = "Maximum Loading on component 1")
plotLoadings(final_model, comp = 2, method = 'mean', contrib = 'max', 
             size.name = 0.8, ndisplay = 30, title = "Maximum Loading on component 2")


# AUROC
auroc(final_model, roc.comp = 1)
auroc(final_model, roc.comp = 2)




#------------------------------------------------------
# 5. Barplots for VIP > 1 Metabolites (Component 1 & 2)
#------------------------------------------------------

process_component <- function(comp_number) {

  # Extract sample coordinates on component 1
  scores <- final_model$variates$X[, comp_number]

  # Metabolite loading direction to group
  loadings <- final_model$loadings$X[, comp_number]

  # Get average Component position for each group
  mean_scores <- data.frame(Group = Y, Score = scores) %>%
    group_by(Group) %>%
    summarise(mean_score = mean(Score)) %>%
    arrange(mean_score)

  # Infer group direction: highest group along comp_number
  # Maximally separated on the component.
  max_group <- as.character(mean_scores$Group[which.max(mean_scores$mean_score)])
  min_group <- as.character(mean_scores$Group[which.min(mean_scores$mean_score)])

  # Convert the factor levels to character strings before using them in ifelse:
  max_group_chr <- as.character(max_group)
  min_group_chr <- as.character(min_group)
  
  metabolite_contrib_group <- ifelse(loadings > 0, max_group_chr, min_group_chr)
  
  metabolite_info <- data.frame(
    Code = names(loadings),
    Loading = loadings,
    Group = factor(metabolite_contrib_group, levels = levels(Y))
  ) %>%
    mutate(IonMode = case_when(
      startsWith(Code, "pos") ~ "Positive",
      startsWith(Code, "neg") ~ "Negative",
      TRUE ~ "Unknown"
    ))

  # Filter by VIP > 1 in comp_number
  vip_scores <- vip(final_model)
  vip_filter <- names(which(vip_scores[, comp_number] > 1))
  vip_df <- metabolite_info %>% filter(Code %in% vip_filter) %>%
  mutate(IonMode = ifelse(startsWith(Code, "pos"), "Positive", 
                          ifelse(startsWith(Code, "neg"), "Negative", "Unknown"))) 
  
  # Create summary
  vip_summary <- vip_df %>%
  group_by(Group, IonMode) %>%
  summarise(n_metabolites = n_distinct(Code), .groups = "drop")

  # Print message and summary
  message("🟢 Number of VIP > 1 metabolites by group and ion mode:")
  print(vip_summary)

  # Bar plot
  plot_obj <- ggplot(vip_df, aes(x = reorder(Code, Loading), y = Loading, fill = Group)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = paste("Component", comp_number, ": Loadings (VIP > 1)"),
      x = "Metabolites", y = "Maximum Loading"
    ) +
    scale_fill_manual(values = group_colors) +
    theme_bw(base_size = 8) +
    theme(axis.text.y = element_text(size = 6))
  
  vip_annotated <- vip_df %>% left_join(metabo_descp, by = c("Code" = "Code"))
  
  return(list(plot = plot_obj, vip_annotated = vip_annotated))
}

# Generate for Component 1 & 2
comp1_out <- process_component(1)
comp2_out <- process_component(2)

# Combine Plots
grid.arrange(
  comp1_out$plot, indv, comp2_out$plot,
  layout_matrix = rbind(c(1, 2), c(1, 3)),
  heights = c(1.5, 2)
)

# Save Annotated Tables
comp1_tuned_VIPcutoff_desp <- comp1_out$vip_annotated
comp2_tuned_VIPcutoff_desp <- comp2_out$vip_annotated

