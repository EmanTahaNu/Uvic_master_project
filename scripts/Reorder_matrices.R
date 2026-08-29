target_cell_lines <- matched_cell_lines$encode_name

print(target_cell_lines)



#`expression matrix`
mat_log_aligned <- readRDS(file.path(MATRICES_DIR, "expression_matrix_aligned_log2TPM.rds"))

missing_exp <- setdiff(target_cell_lines, colnames(mat))

for (cl in missing_exp) {
  mat[, cl]  <- NA
  mat_log[, cl] <- NA
}

#reorder columns

mat_aligned  <- mat[,  c("gene_id", target_cell_lines)]
mat_log_aligned <- mat_log[,  c("gene_id", target_cell_lines)]

# save files

saveRDS(mat_aligned,     file.path(MATRICES_DIR, "expression_matrix_aligned_TPM.rds"))
saveRDS(mat_log_aligned, file.path(MATRICES_DIR, "expression_matrix_aligned_log2TPM.rds"))
write_tsv(mat_aligned,   file.path(MATRICES_DIR, "expression_matrix_alignedTPM.tsv"))
write_tsv(mat_log_aligned,   file.path(MATRICES_DIR, "expression_matrix_aligned_log2TPM.tsv"))
#----------------------------------------

# Histone matrices

for (this_mark in top_marks) {
  h_mat <- readRDS(file.path(MATRICES_DIR,
                             paste0("histone_matrix_", this_mark, ".rds")))
  missing_hist <- setdiff(target_cell_lines, colnames(h_mat))

  for (cl in missing_hist) {
    h_mat[, cl] <- NA
  }
  
# reorder
 h_mat_aligned <- h_mat[, c("cCRE_id", target_cell_lines)]
  
# save as NEW files
 saveRDS(h_mat_aligned, file.path(MATRICES_DIR,
                                   paste0("histone_matrix_", this_mark, "_aligned.rds")))
 write_tsv(h_mat_aligned, file.path(MATRICES_DIR,
                                     paste0("histone_matrix_", this_mark, "_aligned.tsv.gz")))
  
}
  
#-----------------------------------  
 # Metabolites matrix 

# reload fresh metabolite matrix
metabolite_raw <- read_csv("/Users/emantaha/Downloads/CCLE_metabolomics_20190502.csv")

# keep only matched cell lines
metabolite_raw <- metabolite_raw %>%
  filter(CCLE_ID %in% matched_cell_lines$ccle_id)

# reshape rows = metabolites, columns = cell lines (CCLE names)
metabolite_matrix <- metabolite_raw %>%
  dplyr::select(-DepMap_ID) %>%
  tibble::column_to_rownames("CCLE_ID") %>%
  t() %>%
  as.data.frame()

print(colnames(metabolite_matrix))

# rename CCLE names to ENCODE names using crosswalk
name_map <- setNames(matched_cell_lines$encode_name,
                     matched_cell_lines$ccle_id)

colnames(metabolite_matrix) <- name_map[colnames(metabolite_matrix)]

print(colnames(metabolite_matrix))

# check which of the 26 target cell lines are missing
target_cell_lines <- matched_cell_lines$encode_name
missing_met <- setdiff(target_cell_lines, colnames(metabolite_matrix))
print(missing_met)

# add missing as NA columns
for (cl in missing_met) {
  metabolite_matrix[, cl] <- NA
}

# reorder to match target order (same as histone and expression)
met_aligned <- metabolite_matrix[, target_cell_lines]


# confirm same column order as expression matrix
exp_cols <- colnames(mat_aligned)[colnames(mat_aligned) != "gene_id"]

# save
write_csv(met_aligned %>% rownames_to_column("metabolite"),
          file.path(MATRICES_DIR, "metabolite_matrix_aligned.csv"))


  
  
  
  
  