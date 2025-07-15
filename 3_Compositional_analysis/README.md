# 🧬 Compositional Profiles (Taxonomic Abundance)

This script summarizes microbial taxonomic abundances from a `phyloseq` object (`ps`) at multiple taxonomic ranks and exports both:

- ✅ Raw counts  
- ✅ Relative abundances

---

## 📂 Input

- `phyloseq` object (`ps`) loaded from `.RData`
- `otu_table` and `tax_table` must be present in `ps`
- `sample_data` must contain metadata including a `SampleID` column

---

## 🔍 Script Workflow

### 1. Extract ASV Table

- Extracts the ASV count table and taxonomy table from the `phyloseq` object.

### 2. Sort by Total Abundance

- Ranks ASVs by total abundance across all samples in descending order.

### 3. Merge with Taxonomy

- Adds taxonomic labels (Phylum, Class, Order, Family, Genus) to each ASV.

### 4. Group by Taxonomic Rank

- Summarizes total counts at each level using `dplyr::group_by()` and `summarise()`.

### 5. Merge with Metadata

- Combines sample-level metadata (`sample_data(ps)`) with each taxonomic count table.

### 6. Compute Relative Abundances

- Calculates per-sample relative abundances (each row sums to 1).

### 7. Export Outputs

- For each taxonomic level, the script exports:
  - Raw count tables: `*_Counts_*.xlsx`
  - Relative abundance tables: `*_Rel_*.xlsx`

---

## 📤 Output Files

| File Name               | Description                     |
|------------------------|---------------------------------|
| `1_Counts_Phylum.xlsx` | Raw Phylum-level counts         |
| `1_Rel_Phylum.xlsx`    | Phylum-level relative abundances|
| `2_Counts_Class.xlsx`  | Raw Class-level counts          |
| `2_Rel_Class.xlsx`     | Class-level relative abundances |
| `3_Counts_Order.xlsx`  | Raw Order-level counts          |
| `3_Rel_Order.xlsx`     | Order-level relative abundances |
| `4_Counts_Family.xlsx` | Raw Family-level counts         |
| `4_Rel_Family.xlsx`    | Family-level relative abundances|
| `5_Counts_Genus.xlsx`  | Raw Genus-level counts          |
| `5_Rel_Genus.xlsx`     | Genus-level relative abundances |

---

## 📌 Notes

- Any missing taxonomy (e.g., `NA` in Genus) is labeled as `"unclassified"` to ensure compatibility.
- Prefixes such as `p__`, `c__`, `o__`, `f__`, `g__` are removed from column names for cleaner presentation.
- Row sums of relative abundance tables are verified to equal 1 (within floating-point precision).

---

## ✅ Reproducibility Tips

- Always load the correct `.RData` file containing the `phyloseq` object.
- Use `set.seed()` if you perform any downstream random-based operations.
- Ensure `SampleID` is consistent and unique in `sample_data()` to avoid merge issues.

---



# 📊 Differential Abundance Analysis Across BMI Groups

This section analyzes differences in microbial abundance across BMI groups (Normal weight, Overweight, Obese) at various taxonomic ranks using non-parametric statistical tests.

---

## 📂 Input

- Relative abundance tables with metadata (`*_with_metadata_rel`)
- Sample metadata (`sample_metadata`) with `SampleID` and BMI group labeled as `Group` (factor: N, OW, OB)

---


## 🔍 Analysis Workflow

### 🔹 Step 1: Data Preparation

- Extract taxon abundance (excluding `SampleID`, `Group`) for each rank.
- Ensure `Group` is set as a factor with levels: `N`, `OW`, `OB`.

### 🔹 Step 2: Kruskal-Wallis Test

- Non-parametric test (`kruskal.test`) used to detect significant differences in abundance across BMI groups.
- Significant taxa (p < 0.05) are extracted for further testing.
- Results are exported as:

| File Name                | Description                       |
|-------------------------|-----------------------------------|
| `1_krus_phylum_pval.txt`| Kruskal p-values (Phylum level)   |
| `2_krus_class_pval.txt` | Kruskal p-values (Class level)    |
| `3_krus_order_pval.txt` | Kruskal p-values (Order level)    |
| `4_krus_family_pval.txt`| Kruskal p-values (Family level)   |
| `5_krus_genus_pval.txt` | Kruskal p-values (Genus level)    |

### 🔹 Step 3: Dunn’s Post Hoc Test

- For taxa with significant Kruskal results, Dunn’s post hoc tests are performed with BH (Benjamini-Hochberg) correction.
- Significant pairwise comparisons are retained (adjusted p < 0.05).
- Output files:

| File Name                | Description                           |
|-------------------------|---------------------------------------|
| `1_dunn_sig_phylum.txt` | Significant comparisons (Phylum)      |
| `2_dunn_sig_class.txt`  | Significant comparisons (Class)       |
| `3_dunn_sig_order.txt`  | Significant comparisons (Order)       |
| `4_dunn_sig_family.txt` | Significant comparisons (Family)      |
| `5_dunn_sig_genus.txt`  | Significant comparisons (Genus)       |

---

## 📈 Visualization

### 🔹 Boxplots of Significant Taxa (Genus level)

- Boxplots show relative abundance per BMI group for each significant genus.
- Includes individual points (jittered) and faceted layout for each taxon.
- Custom color palette used:
  - **N**: Grey
  - **OW**: Orange (`#FFA500`)
  - **OB**: Dark Red

### 🔹 Heatmap of Significant Genera

- Visualizes significant genera across samples.
- Heatmap rows: genera; columns: samples (colored by Group).
- Spearman distance used for clustering.
- Color scale:
  - Blue → Yellow → Red

---

## 📊 Stacked Barplots

### 🔹 Phylum-Level Barplot

- Relative abundance (mean per group) shown for top 11 phyla.
- Custom palette used for clearer distinction.

### 🔹 Class-Level Barplot

- Displays top 10 most abundant classes + "Others".
- Computed relative abundance row-wise and reshaped for plotting.
- Colors match those in the Phylum plot for consistency.

---

## 📌 Notes

- NA or missing taxa are excluded from boxplots.
- Significant results are based on adjusted p-values.
- Metadata must include `SampleID` and correctly labeled `Group`.

---

## 🧪 Tools Used

- `dplyr`, `FSA`, `ggplot2`, `ggpubr`, `ComplexHeatmap`, `reshape2`
- Kruskal-Wallis and Dunn's Test (non-parametric stats)
- `phyloseq`-generated abundance tables

---
