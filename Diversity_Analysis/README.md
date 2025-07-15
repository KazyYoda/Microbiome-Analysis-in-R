## 📊 Alpha Diversity Analysis

This script estimates and compares alpha diversity indices across groups in your microbiome dataset.

---

### ✅ Methods Used

- **Observed Richness**: Count of observed ASVs
- **Chao1**: Estimated richness accounting for rare taxa
- **Shannon Index**: Evenness + richness
- **Faith's Phylogenetic Diversity (PD)**: Incorporates phylogenetic relationships

---

### 📁 Input

- `Building_Phyloseq.RData`: Phyloseq object containing OTU table, taxonomy, metadata, and tree
- Loaded from `1.Raw_Data/` directory

---

### 📦 Output Files

| File | Description |
|------|-------------|
| `Alpha_diversity.txt` | Table of diversity metrics merged with sample metadata |
| `kruskal_pval.txt` | Kruskal-Wallis test p-values for each metric |
| Plot Output (RStudio) | Faceted boxplots for all diversity indices |

---

### 🧪 Statistical Testing

- **Kruskal-Wallis** test applied to each diversity metric across `Group` levels:  
  `N` (Normal), `OW` (Overweight), `OB` (Obese)
- p-values reported in `kruskal_pval.txt`

---

### 🎨 Visualization

Boxplots show the distribution of each diversity metric across groups with:
- Jittered sample points
- Manual color mapping
- Facets for side-by-side comparison

---

### 📚 Required Packages

Install the following before running:

```r
install.packages(c("car", "picante", "ggplot2", "tidyr", "dplyr", "rio"))
