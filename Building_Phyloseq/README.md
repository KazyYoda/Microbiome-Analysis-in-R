# Microbiome Analysis in R: Building a Phyloseq Object

This script helps you import and structure amplicon sequencing data using the [`phyloseq`](https://joey711.github.io/phyloseq/) package in R. It combines ASV/OTU counts, taxonomy assignments, metadata, and an optional phylogenetic tree into a unified object for downstream ecological and statistical analyses.

---

## 📁 Input Files
Place the following files in the `1.Raw_Data/` directory:

| File | Description |
|------|-------------|
| `feature-table.xlsx` | ASV/OTU table with samples in columns and ASVs in rows. 💡 If starting from feature-table.tsv exported from a BIOM file, open it in Excel, delete the first row (“Constructed from BIOM file”), and save it as feature-table.xlsx.|
| `taxonomy.tsv` | Taxonomic classifications (usually from QIIME 2 or similar). |
| `Metadata.txt` | Sample metadata file (e.g., sample ID, group, treatment). |
| `tree.nwk` | Rooted phylogenetic tree file in Newick format (optional but recommended). |

---

## 🧰 Required Packages

Install dependencies by running the first section of the script. Packages include:

- `phyloseq` (via Bioconductor)
- `biomformat`
- `ape`
- `readxl`
- `dplyr`
- `tibble`
---

## 🚀 Running the Script

1. Open `build_phyloseq.R` in R or RStudio.
2. Adjust the working directory (`setwd(...)`) if needed.
3. Run all sections to load your data.

Once loaded, you’ll have access to all data structures (`otu_table`, `tax_table`, `sample_data`, `phy_tree`) and can construct your `phyloseq` object in the next step.

---

## ✅ Next Steps

Once your data is loaded, proceed to:
- Clean and align data (matching sample names, filtering NA)
- Construct `phyloseq()` object

```r
ps <- phyloseq(otu_table(...), tax_table(...), sample_data(...), phy_tree(...))
