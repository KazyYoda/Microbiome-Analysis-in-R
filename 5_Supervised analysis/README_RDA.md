# 🔬 Microbiome-Metabolite Analysis via Redundancy Analysis (RDA)

This section investigates how gut microbiota at the **genus level** and **metabolite profiles** are associated across different BMI categories (**Normal, Overweight, Obese**), using **Redundancy Analysis (RDA)** in R.

---

## 📌 Workflow Summary

### Step 1: Data Preparation
- Load genus-level counts and filter to 55 metabolomics-matched samples.
- Add pseudocount (+1) to avoid `log(0)`, then apply CLR transformation using `compositions::clr()`.
- Load pre-scaled metabolomics data (`metabo_scaled`) from `sPLSDA.RData`.

### Step 2: Redundancy Analysis (RDA)
- Fit the RDA model using:
  - **Response:** `metabo_scaled` (log-transformed and scaled)
  - **Explanatory:** `Group` from metadata (BMI group)
- Evaluate:
  - Model summary and adjusted R²
  - Canonical coefficients
  - Variance partitioning: constrained vs. unconstrained

### Step 3: Statistical Testing
- Run permutation-based ANOVA using `vegan::anova.cca()`:
  - Overall model significance
  - Individual constrained axes (e.g., RDA1, RDA2)
  - Group term significance

### Step 4: Visualization
- Use a custom `rda_plot()` function to:
  - Plot samples (colored by BMI group)
  - Plot top contributing features (metabolites/genus)
  - Color arrows by feature type:
    - Genus: **purple** (`#660099`)
    - Positive ions: **green** (`#228B22`)
    - Negative ions: **red** (`#e02b35`)
  - Supports:
    - **Scaling 1:** Sample distances
    - **Scaling 2:** Feature correlations
- Combine plots using `gridExtra::grid.arrange()`

---

## 📊 Statistical Interpretation

- RDA quantifies how much variance in the multi-omics dataset is explained by BMI group.
- Permutation-based ANOVA tests:
  - Overall model significance
  - Significance of RDA axes (e.g., RDA1, RDA2)
  - Group effect (BMI categories: **N**, **OW**, **OB**)
- Adjusted R² provides a corrected measure of explanatory power and prevents overfitting.

---

## 📈 Visualization Outputs

- `rda_plot()` produces informative biplots with:
  - Sample points (colored by BMI group)
  - Top features as arrows (arrow length reflects contribution)
  - Switchable between:
    - **Scaling 1** — emphasizes sample distances
    - **Scaling 2** — emphasizes feature correlations

**Arrow Colors**
| Feature Type   | Color Code |
|----------------|------------|
| Genus          | `#660099`  |
| Positive Ions  | `#228B22`  |
| Negative Ions  | `#e02b35`  |

---

## 📤 Output Files

| Output Name                | Description                                           |
|----------------------------|-------------------------------------------------------|
| `rda_model`                | RDA model object containing fitted multivariate model |
| `rda_plot(scaling = 2)`    | Biplot emphasizing relationships between features     |
| `rda_plot(scaling = 1)`    | Biplot emphasizing sample distances                   |
| `anova_results.txt` (opt.) | Significance results from permutation ANOVA          |
| `RDA1_vs_RDA2_biplot.pdf`  | Exported ggplot2 visualization of RDA                |
