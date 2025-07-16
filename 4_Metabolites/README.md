# Microbiome Metabolites Analysis

This project analyzes the relative abundance and regulation of positive- and negative-ion metabolites across three groups: Normal (N), Overweight (OW), and Obese (OB). It uses log2 and log2 fold-change calculations, summary statistics, and visualizations including heatmaps and bar plots.

---

## 📁 Directory Structure

```
4.Metabolites/
├── POS_Metabo.xlsx          # Positive-ion metabolite intensity data
├── NEG_Metabo.xlsx          # Negative-ion metabolite intensity data
├── Metabolite_analysis.R      # Main R script for analysis and visualization
├── output/
│   ├── 1_log2FC_pos_reltoN.xlsx
│   ├── 2_log2FC_neg_reltoN.xlsx
│   ├── 1_log2FC_grouped_pos.xlsx
│   └── 2_log2FC_grouped_neg.xlsx
```

---

## 🛠️ Requirements

Install the required R packages:

```r
install.packages(c("dplyr", "tidyr", "readxl", "ggplot2", "car",
                   "ComplexHeatmap", "circlize"))
```

---

## 🔄 Workflow Summary

### Step 1: Data Import and Metadata Preparation

* Import metabolite intensity tables (POS and NEG ion modes).
* Filter and format metadata to match 55 relevant samples.

### Step 2: Log2 Fold Change Calculation

* Calculate mean metabolite intensity in Normal group (baseline).
* Compute log2 fold-change (log2FC) relative to baseline using `sweep()`.
* Save output tables for both POS and NEG datasets.

### Step 3: Summary Statistics

* Summarize log2FC per group per metabolite.
* Count number of upregulated/downregulated metabolites (relative to N).
* Export grouped statistics.

### Step 4: Visualization

* Generate heatmaps (Clustered by Spearman correlation).
* Plot grouped bar charts with error bars (SD) showing mean log2FC.

---

## 📊 Output Files

* **log2FC\_pos\_reltoN.xlsx**: Log2FC values for each positive-ion metabolite.
* **log2FC\_grouped\_pos.xlsx**: Mean log2FC and SD by group (POS).
* **log2FC\_neg\_reltoN.xlsx**: Log2FC values for each negative-ion metabolite.
* **log2FC\_grouped\_neg.xlsx**: Mean log2FC and SD by group (NEG).

---

## 📎 Notes

* Be sure your `POS_Metabo.xlsx` and `NEG_Metabo.xlsx` files contain matching `SampleID` columns.
* Heatmaps are based on log2 values and colored according to relative abundance.
* Metabolites are renamed as `pos1`, `pos2`, ... or `neg1`, `neg2`, ... to maintain consistency.

---

