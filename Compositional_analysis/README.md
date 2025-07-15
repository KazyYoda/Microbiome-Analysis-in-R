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
