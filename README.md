# Mantle Cell Lymphoma Affymetrix Microarray Pipeline

This repository contains R scripts for preprocessing Affymetrix 3′ expression microarrays from **mantle cell lymphoma (MCL)** samples. The workflow performs RMA normalization, maps probe IDs to gene symbols, and merges multiple batches into a single gene-level expression matrix suitable for downstream biomarker analyses (e.g., ROR1 expression ranking across MCL cases).

---

## Overview of the Workflow

1. **Batch-wise preprocessing of Affymetrix `.CEL` files**
   - Loops over user-specified subfolders (e.g., `1`, `3`, `5`, `7`, `8`) containing raw MCL `.CEL` files.
   - Uses `affy::ReadAffy()` and `rma()` to generate log2, RMA-normalized expression values.
   - Annotates probe sets to gene symbols via `hgu133plus2.db`.
   - Exports one tab-delimited file per batch with columns:
     - `ProbeID`
     - `GeneSymbol`
     - One column per MCL sample.

2. **Merging batches by GeneSymbol**
   - Reads all per-batch text files into R.
   - Merges them on `GeneSymbol` using `dplyr::full_join()` (configurable to `inner_join()` if only shared genes are needed).
   - Produces a single, combined MCL expression matrix (`MCL_microarray_combined_by_gene_symbol.txt`).

3. **Downstream analysis (example: ROR1 in MCL)**
   - The combined matrix can be used to extract **ROR1** and other genes of interest.
   - Example downstream scripts (not shown here) can generate barplots of ROR1 expression across MCL samples, rank-ordered from lowest to highest, as in the illustrative figure.

---

## Directory Structure

Example layout used in the scripts:

```text
/Users/yourname/Downloads/MCL/WorkDir
├── 1/                      # Raw MCL .CEL files (batch 1)
├── 3/                      # Raw MCL .CEL files (batch 3)
├── 5/
├── 7/
├── 8/
└── Expression files/       # RMA- and symbol-annotated outputs
    ├── MCL_microarray_gene_symbols_1.txt
    ├── MCL_microarray_gene_symbols_3.txt
    ├── MCL_microarray_gene_symbols_5.txt
    ├── MCL_microarray_gene_symbols_7.txt
    └── MCL_microarray_gene_symbols_8.txt
```

You can adjust `base_dir`, `subfolders`, and the paths to match your own MCL dataset organization.

---

## Dependencies

The scripts rely on the following R packages:

- **Bioconductor**
  - `affy`
  - `annotate`
  - `hgu133plus2.db` (for GPL570 / HG-U133 Plus 2.0 MCL arrays)
- **CRAN**
  - `dplyr`
  - `purrr`

Installation is handled automatically in the script using `BiocManager::install()` and `install.packages()` if a package is missing.

If you are using newer ST-type Affymetrix arrays, you can adapt the commented `oligo`-based code and the appropriate `pd.*` annotation package, but the current repository is tuned to classic 3′ expression arrays used in MCL datasets.

---

## Usage

### 1. Preprocess MCL `.CEL` files

In the first section of the script:

```r
base_dir  <- "/Users/yourname/Downloads/MCL/WorkDir"
subfolders <- c("1", "3", "5", "7", "8")
```

- Set `base_dir` to the root directory containing your MCL `.CEL` subfolders.
- Set `subfolders` to the names of the folders you want to process.
- Run the loop. For each folder, the script:
  - Reads all `.CEL` files.
  - Applies RMA normalization.
  - Maps probes → gene symbols.
  - Writes `MCL_microarray_gene_symbols_<subfolder>.txt` into the chosen output directory.

### 2. Combine batch-level files by GeneSymbol

In the second section:

```r
setwd("/Users/yourname/Downloads/MCL/WorkDir/Expression files")
files <- c(
  "MCL_microarray_gene_symbols_1.txt",
  "MCL_microarray_gene_symbols_3.txt",
  "MCL_microarray_gene_symbols_5.txt",
  "MCL_microarray_gene_symbols_7.txt",
  "MCL_microarray_gene_symbols_8.txt"
)
```

Then:

```r
list_dfs <- lapply(files, function(f) {
  read.delim(f, header = TRUE, sep = "	", stringsAsFactors = FALSE)
})

combined_df <- purrr::reduce(
  list_dfs,
  function(x, y) dplyr::full_join(x, y, by = "GeneSymbol")
)

write.table(
  combined_df,
  file = "MCL_microarray_combined_by_gene_symbol.txt",
  sep = "	",
  quote = FALSE,
  row.names = FALSE
)
```

This will create a single gene-level matrix for all MCL samples across batches.

### 3. Downstream analyses

With `MCL_microarray_combined_by_gene_symbol.txt`:

- Load the matrix into R or Python.
- Subset rows where `GeneSymbol == "ROR1"` or any other biomarker panel.
- Generate MCL-focused visualizations (e.g., ROR1 expression distribution, heatmaps, clustering, survival modeling).

---

## Notes

- Rows with missing `GeneSymbol` can be filtered out if desired before merging to simplify downstream analyses.
- Sample column names can be cleaned (e.g., removing `.CEL` suffixes) for more readable plots.
- The current setup is intentionally specific to **mantle cell lymphoma** microarray datasets; extension to other diseases simply requires adjusting file paths and metadata but does not change the core normalization/merge logic.

---

## Contact

For questions or suggestions related to the MCL pipeline:

- **Maintainer:** Fahad S. Ahmed, MD  
- **Email:** fahadshabbirahmed@gmail.com
