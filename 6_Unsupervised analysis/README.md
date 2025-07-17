# Microbiome-Metabolite Multi-Omics Integration using Multiple Factor Analysis (MFA)

This script performs integrative analysis of gut microbiota (genus-level) and metabolomics data (positive and negative ions) using Multiple Factor Analysis (MFA). The goal is to uncover shared structures and group differences (e.g., BMI categories) across datasets in a reduced dimensional space.

---

## 📂 Input Files

| File Name                                  | Description                                  |
|--------------------------------------------|----------------------------------------------|
| `5_Counts_Genus.xlsx`                      | Genus-level count data                       |
| `5_krus_genus_pval.txt`                    | Kruskal-Wallis p-values for genus filtering  |
| `sPLSDA.RData`                             | Contains `metabo_scaled` and `Metadata_55`  |

---

## 📦 Required R Packages

```r
library(dplyr)
library(readxl)
library(compositions)
library(FactoMineR)
library(factoextra)
library(ggplot2)
```


# ⚙️ Workflow Summary

## Step 1: Data Preparation

- Load genus-level counts and filter to **55 metabolomics-matched samples**.
- Add **pseudocount (+1)** to avoid log(0), then apply **CLR transformation** using `compositions::clr()`.
- Load **pre-scaled metabolomics data (`metabo_scaled`)** from `sPLSDA.RData`.

## Step 2: Redundancy Analysis (RDA)

- **Fit the RDA model** using:
  - **Response**: `metabo_scaled` (log-transformed + scaled)
  - **Explanatory**: `Group` from metadata (BMI group)
- **Evaluate**:
  - Model summary and **adjusted R²**
  - **Canonical coefficients**
  - **Variance partitioning**: constrained vs. unconstrained

## Step 3: Statistical Testing

- Run **permutation-based ANOVA** using `vegan::anova.cca()`:
  - Overall model
  - Individual constrained axes (e.g., RDA1, RDA2)
  - Group term significance

## Step 4: Visualization

- Use custom `rda_plot()` function to:
  - Plot **samples by BMI group**
  - Plot **top contributing features** (metabolites/genus)
  - **Arrow coloring by feature type**: `Genus`, `positive_ion`, `negative_ion`
  - Supports:
    - **Scaling 1**: Emphasizes sample distances
    - **Scaling 2**: Emphasizes feature correlations
  - Combine plots with `gridExtra::grid.arrange()`

---

# 📊 Statistical Interpretation

- RDA quantifies how much variance in the **multi-omics dataset** is explained by the **BMI group**.
- **Permutation-based ANOVA** tests:
  - Overall model significance
  - Significance of RDA axes (e.g., RDA1, RDA2)
  - **Group effect** (BMI categories: N, OW, OB)
- **Adjusted R²** prevents overfitting and reflects true explanatory power.

---

# 📈 Visualization Outputs

### `rda_plot()` generates:

- Sample points **colored by BMI group**
- **Top features as arrows** (length reflects contribution)
- Switchable Scaling:
  - **Scaling 1**: Sample distances
  - **Scaling 2**: Feature correlations

### Arrow Color by Feature Type:

| Feature        | Color       | Hex       |
|----------------|-------------|-----------|
| Genus          | Purple      | `#660099` |
| Positive Ions  | Green       | `#228B22` |
| Negative Ions  | Red         | `#e02b35` |

---

# 📤 Output Files

| Output Name               | Description                                                |
|---------------------------|------------------------------------------------------------|
| `rda_model`               | RDA model object containing fitted multivariate regression |
| `rda_plot(scaling = 2)`   | Biplot emphasizing feature correlations                    |
| `rda_plot(scaling = 1)`   | Biplot emphasizing sample distances                        |
| `anova_results.txt` (opt) | Permutation-based ANOVA significance test results          |
| `RDA1_vs_RDA2_biplot.png` | Exported `ggplot2` visualization of the RDA analysis       |
