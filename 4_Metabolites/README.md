# Microbiome Metabolites Analysis

This script analyzes metabolite profiles (positive and negative ion modes) in relation to obesity status across three groups: Normal (N), Overweight (OW), and Obese (OB). It computes log2 and log2 fold changes, summarizes results, and visualizes differential patterns using heatmaps and bar plots.

---

## 📥 Input

* `POS_Metabo.xlsx` — Positive-ion metabolite data
* `NEG_Metabo.xlsx` — Negative-ion metabolite data

Each file must include:

* `SampleID` column
* `Group` column (values: N, OW, OB)
* Metabolite columns 

---

## 📦 Required Packages

Install the required R packages:

```r
install.packages(c("dplyr", "tidyr", "readxl", "ggplot2", "car",
                   "ComplexHeatmap", "circlize"))
```

---

## 🔄 Workflow Summary

### Step 1: Data Import and Metadata Preparation

* Load both POS and NEG metabolite datasets.
* Ensure `SampleID` alignment between POS and NEG.
* Extract `SampleID` and `Group` metadata for 55 samples.

### Step 2: Compute Log2 Fold Change (log2FC)

* Calculate mean metabolite abundance in Normal (N) group.
* Apply `sweep()` to compute log2FC relative to Normal group.
* If data contain `NA` or zero values then add 1 to numeric columns only to avoid log2(0)
```r
data <- data %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0) + 1))
```

### Step 3: Data Summary 

* Convert data to long format.
* Calculate mean ± SD log2FC for each metabolite within each group.
* Classify metabolites as Upregulated / Downregulated / No Change.
* Count total regulated metabolites per group.

### Step 4: Visualization

* Generate clustered heatmaps (Spearman distance):
  * log2-transformed metabolite realtive abundance
    
* Create grouped bar plots with error bars (SD), colored by direction.
  * log2FC

### Step 5: Export Outputs

* Export log2FC tables and grouped summaries as `.xlsx` files for further use.

---

## 📊 Statistical Testing and Interpretation

* Interpretation is based on magnitude and direction of log2FC:

  * Positive = Upregulated (log2FC > 0)
  * Negative = Downregulated (log2FC < 0)
  * No change = ≈ 0

* The analysis includes non-parametric and parametric tests to assess metabolite
  differences between groups:
> **Note**: Normality and homogeneity of variance should be assessed prior to performing statistical analysis.

### 1. **Kruskal-Wallis Test**
- A non-parametric method used to detect differences in median values across more than two groups.
- Suitable when the data is not normally distributed or variances are unequal.
- **Null hypothesis (H₀)**: All groups come from the same distribution.
- **Decision Rule**: Reject H₀ if `p < 0.05`.

### 2. **Dunn’s Post Hoc Test**
- Conducted after a significant Kruskal-Wallis result.
- Performs pairwise comparisons between groups.
- Adjusted p-values are applied using methods like Bonferroni or Benjamini-Hochberg correction.

### 3. **One-Way ANOVA**
- A parametric method to compare mean differences across groups.
- Assumes normally distributed data and equal variances (homoscedasticity).
- **Null hypothesis (H₀)**: All group means are equal.
- **Decision Rule**: Reject H₀ if `p < 0.05`.

### 4. **Tukey’s HSD Test**
- Post hoc analysis following a significant ANOVA result.
- Identifies which specific groups differ from each other.

---

## 📈 Visualization

* **Heatmaps**:

  * Clustered by Spearman correlation
  * Colored using diverging palette centered at 0 (gray for no change)

* **Bar Plots**:

  * Mean log2FC ± SD
  * Colored by regulation direction
  * Group-specific facet panels

---

## 📁 Output Files

| File Name                   | Description                                        |
| --------------------------- | -------------------------------------------------- |
| `1_log2FC_pos_reltoN.xlsx`  | Log2FC of POS metabolites relative to Normal group |
| `2_log2FC_grouped_pos.xlsx` | Mean ± SD log2FC of POS metabolites per group      |
| `2_log2FC_neg_reltoN.xlsx`  | Log2FC of NEG metabolites relative to Normal group |
| `2_log2FC_grouped_neg.xlsx` | Mean ± SD log2FC of NEG metabolites per group      |

---

## 📎 Notes

* Make sure both datasets (`POS_Metabo.xlsx` and `NEG_Metabo.xlsx`) share identical `SampleID` ordering.
* `log2FC` values are calculated relative to **Normal (N)** group.
* Zero or missing values should be handled prior to transformation.
* Metabolites are renamed to `pos1`, `pos2`, ..., `neg1`, `neg2`, ... for plotting clarity.

---
