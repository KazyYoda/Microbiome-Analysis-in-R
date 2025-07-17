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



# ⚙️ Workflow Summary: Multiple Factor Analysis (MFA)

## Step 1: Data Preparation
- Load genus count table and filter genera with p < 0.05 from Kruskal-Wallis test.
- Apply CLR transformation on filtered genus counts (with pseudocount +1).
- Load pre-scaled metabolomics data (`metabo_scaled`) from prior sPLS-DA step.
- Merge CLR-transformed genus data with metabolomics (positive + negative ion) and group metadata.

## Step 2: Perform MFA
- Use `FactoMineR::MFA()` with group structure:
  - Group: BMI Group (qualitative)
  - Positive Ions (95 features)
  - Negative Ions (38 features)
  - Genus (26 features)
- Configure group types (`n`, `c`, `c`, `c`) and exclude “Group” from active analysis (`num.group.sup = 1`).
- Generate scree plot for eigenvalue inspection.
- type: the type of variables in each group. By default, all variables are quantitative and scaled to unit variance. Allowed values include:
“c” or “s” for quantitative variables. If “s”, the variables are scaled to unit variance.
“n” for categorical variables.
“f” for frequencies (from a contingency tables).
- Reference: MFA - Multiple Factor Analysis in R (https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/116-mfa-multiple-factor-analysis-in-r-essentials/)

## Step 3: Group Variable Contributions
- Use `get_mfa_var()` to extract:
  - Coordinates
  - Cos² (representation quality)
  - Contributions to each dimension
- Visualize:
  - Group contributions to Dimensions 1–2 and 3–4
  - Barplots for group contributions to specific axes

## Step 4: Sample Projections (Individuals)
- Use `fviz_mfa_ind()` to visualize sample projections colored by BMI group (N, OW, OB).
- Add confidence ellipses per group.
- Customize aesthetic: axis labels, legend, and spacing.

## Step 5: Quantitative Variable Analysis
- Plot top 50 contributing variables across Dim 1–4 using `fviz_contrib()`.
- Display correlation maps for variables on Dimensions 1–2 and 3–4.
- Arrows reflect contribution strength and direction in feature space.

## Step 6: Dimension Descriptions
- Use `dimdesc()` to describe the most contributing variables (features) to each dimension.
- Helps interpret biological drivers of MFA axes.

---

## 📊 Statistical Interpretation
- MFA integrates multiple datasets while preserving their structure and balance.
- Dimension 1 & 2 capture the most shared variance across microbiota and metabolomics.
- Group effects (e.g., BMI) are observed in projections.
- Variable contributions highlight which genera/metabolites drive these differences.
- Cos² and contribution scores guide selection of key features for follow-up.

---

## 📈 Visualization Outputs

| Plot Type                   | Description                                                   |
|----------------------------|---------------------------------------------------------------|
| Scree Plot                 | Eigenvalues showing variance explained by each dimension      |
| Group Contribution Barplots| Contributions of genus, pos/neg ions to each axis             |
| Group Correlation Maps     | Variable group locations on MFA dimensions                    |
| Individual Sample Biplots  | Sample clustering colored by group with ellipses              |
| Variable Contribution      | Feature-wise contribution barplots for dimensions 1–4         |
| Variable Correlation Maps  | Quantitative feature correlations shown as vectors            |

---

## 📤 Output Files

| Output Name           | Description                                                     |
|-----------------------|-----------------------------------------------------------------|
| `res.mfa`             | MFA result object containing scores and loadings                |
| `fviz_mfa_ind()` plot | Sample projection on MFA space with group ellipses              |
| `fviz_mfa_var()` plots| Variable contributions and group correlations                   |
| `fviz_contrib()` plots| Contributing features per dimension                         |
| `dimdesc(res.mfa, axis)`| Lists significant variables for each MFA axis                  |

---

## 📝 Notes
- Pseudocounts (+1) are used before CLR transformation to handle zeros.
- Sample order is validated (`identical()`) to prevent alignment errors.
- Feature colors can be customized.
- Be sure to load `sPLSDA.RData` to access `metabo_scaled` and `Metadata_55`.
