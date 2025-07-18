# Microbiome Analysis: Metabolite Classification with sPLS-DA

This R script performs a **sparse Partial Least Squares Discriminant Analysis (sPLS-DA)** on positive and negative ion metabolite datasets obtained from microbiome samples. It identifies key discriminative metabolites between defined groups (e.g., N, OW, OB), visualizes sample separation, and extracts variable importance metrics.

---

## 📥 Input

| File Name                        | Description                                      |
|----------------------------------|--------------------------------------------------|
| `posmetabo`, `negmetabo`        | Normalized positive and negative ion metabolite data (loaded from workspace) |
| `Metadata_55`                   | Sample metadata including group assignments      |
| `metabo_descp.xlsx`             | Mapping table: metabolite codes → names/classes  |

---

## 📦 Required Packages

```r
library(mixOmics)
library(MASS)
library(lattice)
library(ggplot2)
library(dplyr)
library(readxl)
library(gridExtra)
```

# ⚙️ sPLS-DA Metabolomics Workflow

This script outlines a complete sPLS-DA pipeline in R for multivariate analysis of metabolomics data, including model optimization, visualization, and statistical interpretation.

---

## 🚀 Script Workflow

### 1. Data Preparation
- Combine `posmetabo` and `negmetabo` using `cbind()`
- Log2-transform and scale the metabolite matrix
- Extract group labels from `metadata`

### 2. Basic sPLS-DA
- Run initial `splsda()` with `ncomp = 2`
- Generate individual and variable plots
- Retrieve VIP scores

### 3. Model Tuning
- Use cross-validation to:
  - Choose optimal number of components
  - Tune `keepX` (number of features to retain per component)

### 4. Final Model & Visualization
- Build final model using optimal parameters
- Generate plots:
  - sPLS-DA score plots with ellipses
  - VIP > 1 barplots for Component 1 and 2
  - Loading plots for top-contributing metabolites
  - AUROC curve plots per component

### 5. Interpretation & Group Contribution
- Match metabolite loading direction to sample group separation
- Identify key metabolites with VIP > 1
- Group-wise and ion-mode-wise summaries
- Annotate results using `metabo_descp`

---

## 📊 Statistical Testing and Model Evaluation

- **Cross-Validation**:  
  5-fold, 100 repeats to evaluate Balanced Error Rate (BER) and AUROC.
  
- **Model Optimization**:
  - `perf()` used to assess component performance
  - `tune.splsda()` selects optimal `keepX`
  
- **AUROC**:  
  Computed per component to quantify classification performance.
  
- **VIP Scores**:  
  Threshold of VIP > 1 used to define high-contributing features.

> **Note**: Normality and homogeneity of variance should be assessed prior to performing statistical analysis.

---

## 🖼️ Visualizations

- **Score Plots**: Colored ellipses and group separation via `ggplot2`
- **Loading Barplots**: VIP > 1 metabolites with group contribution
- **AUROC Curves**: Evaluate model performance

---

## 📁 Output Files

| File / Object Name              | Description                                               |
|--------------------------------|-----------------------------------------------------------|
| `final_model`                  | Final tuned sPLS-DA model object                          |
| `vip_scores`                   | VIP score matrix by component                            |
| `indv`                         | ggplot2-based score plot                                 |
| `met_comp1`, `met_comp2`       | Barplots of Component 1 and 2 metabolite loadings (VIP > 1) |
| `comp1_tuned_VIPcutoff_desp`   | VIP > 1 metabolites in Component 1 with annotations       |
| `comp2_tuned_VIPcutoff_desp`   | VIP > 1 metabolites in Component 2 with annotations       |

---

## 📝 Notes

- Ensure that zero values are handled before log-transformation (e.g., add pseudocount).
- Metadata group labels (`Metadata_55$Group`) must be factor-type.
- Adjust visual themes and component labels as needed based on the explained variance.
