###########################################################################
##########               Folder 1, 3, 5, 7, 8                 #############
###########################################################################

###############################################################################
# R Script: Loop Over Multiple Folders to Process Affymetrix .CEL Files
###############################################################################

# 1) If not already installed, install needed Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_pkgs <- c("affy", "annotate", "hgu133plus2.db")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, update = FALSE)
  }
}

# 2) Load required libraries
library(affy)
library(annotate)
library(hgu133plus2.db)  # Adjust if your platform requires a different annotation library

# 3) Define the base directory and subfolders to process
base_dir  <- "/Users/User/Downloads/MCL/WorkDir"
subfolders <- c("1", "3", "5", "7", "8")  # Extend or reduce as needed

# 4) Loop through each subfolder
for (subdir in subfolders) {
  
  # Construct the full path to the subfolder
  path_to_cels <- file.path(base_dir, subdir)
  
  # Change the working directory to that folder
  setwd(path_to_cels)
  message("Processing folder: ", path_to_cels)
  
  # 4a) Read .CEL files in this folder
  affy_data <- ReadAffy()  # For classic 3’ arrays, returns an AffyBatch
  
  # 4b) Perform RMA normalization
  eset <- rma(affy_data)   # Produces an ExpressionSet (log2, RMA-normalized)
  
  # 4c) Extract the expression matrix
  expr_matrix <- exprs(eset)
  
  # 4d) Map probe IDs to gene symbols
  probe_ids    <- rownames(expr_matrix)
  gene_symbols <- getSYMBOL(probe_ids, "hgu133plus2") 
  # Adjust key if using a different annotation library
  
  # 4e) Combine into a data frame
  out_df <- data.frame(
    ProbeID    = probe_ids,
    GeneSymbol = gene_symbols,
    expr_matrix,
    check.names = FALSE
  )
  
  # 4f) Write a tab‐delimited file named uniquely for each folder
  out_filename <- paste0("MCL_microarray_gene_symbols_", subdir, ".txt")
  write.table(
    out_df,
    file = out_filename,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  message("  -> Finished. Output file: ", out_filename, "\n")
}

###############################################################################
# Note: If your arrays are modern ST-type arrays (not classic 3' expression arrays),
# use the oligo package. For example:
#
# library(oligo)
# for (subdir in subfolders) {
#   path_to_cels <- file.path(base_dir, subdir)
#   setwd(path_to_cels)
#   celFiles <- list.celfiles(full.names=TRUE)
#   affy_data <- read.celfiles(celFiles)   # Returns ExpressionFeatureSet
#   eset <- rma(affy_data)
#   # Then proceed with annotation & writing out similarly, 
#   # but using a 'pd.*' annotation package if needed.
# }
###############################################################################




###########################################################################
##########                    Combine files.                  #############
###########################################################################

###############################################################################
# R Script: Combine Multiple Microarray Expression Files by GeneSymbol
###############################################################################

# 1) Install required packages if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("purrr", quietly = TRUE)) {
  install.packages("purrr")
}

# 2) Load libraries
library(dplyr)  # for joins
library(purrr)  # for reduce()

# 3) Define the file paths
#    Adjust these to match the actual names/locations of your data files.
setwd("/Users/User/Downloads/MCL/WorkDir/Expression files")
files <- c(
  "MCL_microarray_gene_symbols_1.txt",
  "MCL_microarray_gene_symbols_3.txt",
  "MCL_microarray_gene_symbols_5.txt",
  "MCL_microarray_gene_symbols_7.txt",
  "MCL_microarray_gene_symbols_8.txt"
  # Add or remove filenames as needed, e.g. "MCL_microarray_gene_symbols_3.txt"
)

# 4) Read each file into a list of data frames
list_dfs <- lapply(files, function(f) {
  df <- read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Optional: remove rows with missing GeneSymbol
  # df <- df[!is.na(df$GeneSymbol) & df$GeneSymbol != "", ]
  
  # Optional: rename columns if needed to avoid duplicated sample names, e.g.:
  # colnames(df) <- sub("\\.CEL$", "", colnames(df))
  
  df
})

# 5) Merge all data frames by "GeneSymbol"
#    "full_join()" keeps all GeneSymbols present in any file.
#    If you only want GeneSymbols shared by all files, use "inner_join()" instead.
combined_df <- reduce(list_dfs, function(x, y) {
  full_join(x, y, by = "GeneSymbol")
})

# 6) Write the combined table to a tab-delimited file
#write.table(
#  combined_df,
#  file = "MCL_microarray_combined_by_gene_symbol.txt",
#  sep = "\t",
#  quote = FALSE,
#  row.names = FALSE
#)

#cat("Done! Created 'MCL_microarray_combined_by_gene_symbol.txt'.\n")


