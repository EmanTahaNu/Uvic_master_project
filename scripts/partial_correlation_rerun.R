library(dplyr)
library(readr)
library(ppcor)

CORR_SOURCE_DIR <- file.path(RESULTS_DIR, "pls_correlation_full_rerun")
OUT_DIR         <- file.path(RESULTS_DIR, "partial_correlation_v29")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

top_marks <- c("H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3")

FDR_CORR_THRESHOLD <- 0.1  
MIN_N   <- 4      # minimum complete cell lines 

# Load V29 cCRE-gene assignment

PLS_annotated <- read.delim(
  file.path(RESULTS_DIR, "gene_assignment_v29", "PLS_gene_assignment_v29_pc.txt"),
  stringsAsFactors = FALSE
)

PLS_annotated$gene_id <- sub("\\..*$", "", PLS_annotated$gene_id)   # strip version suffix
PLS_annotated <- PLS_annotated %>%
  dplyr::select(cCRE_id, gene_id) %>%
  distinct()


# Load metabolite & expression matrices

metabolite_matrix <- read.csv(
  file.path(MATRICES_DIR, "metabolite_matrix_aligned.csv"),
  stringsAsFactors = FALSE
)


gene_matrix <- readRDS(file.path(MATRICES_DIR, "expression_matrix_aligned_log2TPM.rds"))
gene_matrix$gene_id <- sub("\\..*$", "", gene_matrix$gene_id)
gene_matrix <- gene_matrix %>% distinct(gene_id, .keep_all = TRUE)


cat("Loaded ", nrow(PLS_annotated), " unique cCRE-gene assignments\n", sep = "")
cat("Loaded ", nrow(metabolite_matrix), " metabolites\n", sep = "")
cat("Loaded ", nrow(gene_matrix), " genes\n", sep = "")

multi_gene_ccres <- PLS_annotated %>% count(cCRE_id) %>% filter(n > 1) %>% nrow()
cat("cCREs with more than one gene assignment: ", multi_gene_ccres, "\n", sep = "")

col_has_data <- function(mat, s) any(!is.na(mat[[s]]))

# STEP 1  —  partial correlation for significant pairs

for (this_mark in top_marks) {
  
  # significant cCRE-metabolite pairs from ordinary correlation
  pls_corr <- read.csv(
    file.path(CORR_SOURCE_DIR, paste0("PLS_correlation_", this_mark, ".csv")),
    stringsAsFactors = FALSE
  )
  sig_pairs <- pls_corr %>% dplyr::filter(spearman_fdr < FDR_CORR_THRESHOLD)
  cat(this_mark, " - correlation-significant pairs (FDR < ",
      FDR_CORR_THRESHOLD, "): ", nrow(sig_pairs), "\n", sep = "")
  
  # histone matrix
  mark_matrix <- readRDS(
    file.path(ALL_MATRICES_DIR, paste0("histone_matrix_", this_mark, "_aligned.rds"))
  )
  
  #  common cell lines: REAL data in all three matrices
  histone_samples    <- setdiff(colnames(mark_matrix), c("cCRE_id", "category"))
  metabolite_samples <- setdiff(colnames(metabolite_matrix), "metabolite")
  expression_samples <- setdiff(colnames(gene_matrix), "gene_id")
  
  named_common <- Reduce(intersect,
                         list(histone_samples, metabolite_samples, expression_samples))
  
  common_samples <- named_common[vapply(named_common, function(s)
    col_has_data(mark_matrix, s) &&
      col_has_data(metabolite_matrix, s) &&
      col_has_data(gene_matrix, s),
    logical(1))]
  
  cat(this_mark, " - named common: ", length(named_common),
      " | with real data: ", length(common_samples), "\n", sep = "")
  cat(this_mark, " - cell lines: ", paste(sort(common_samples), collapse = ", "), "\n", sep = "")
  
  # triplets: significant pairs joined to assigned gene(s)
  test_pairs <- sig_pairs %>%
    dplyr::inner_join(PLS_annotated, by = "cCRE_id") %>%
    distinct(cCRE_id, gene_id, metabolite) %>%
    dplyr::filter(cCRE_id %in% mark_matrix$cCRE_id,
                  gene_id %in% gene_matrix$gene_id)
  
  cat(this_mark, " - cCRE-metabolite-gene triplets: ", nrow(test_pairs),
      " (", dplyr::n_distinct(test_pairs$cCRE_id), " unique cCREs)\n", sep = "")
  
  results <- vector("list", nrow(test_pairs))
  skipped <- vector("list", nrow(test_pairs))
  result_counter <- 0
  skip_counter   <- 0
  
  for (i in seq_len(nrow(test_pairs))) {
    
    this_ccre  <- test_pairs$cCRE_id[i]
    metabolite <- test_pairs$metabolite[i]
    gene       <- test_pairs$gene_id[i]
    
    y <- mark_matrix[mark_matrix$cCRE_id == this_ccre, common_samples, drop = TRUE]
    x <- metabolite_matrix[metabolite_matrix$metabolite == metabolite, common_samples, drop = TRUE]
    z <- gene_matrix[gene_matrix$gene_id == gene, common_samples, drop = TRUE]
    
    if (length(y) != length(common_samples) ||
        length(x) != length(common_samples) ||
        length(z) != length(common_samples)) {
      skip_counter <- skip_counter + 1
      skipped[[skip_counter]] <- data.frame(
        cCRE_id = this_ccre, gene_id = gene, metabolite = metabolite,
        n = NA_integer_, reason = "missing_or_duplicate_matrix_row")
      next
    }
    
    keep <- complete.cases(x, y, z)
    x <- as.numeric(x[keep]); y <- as.numeric(y[keep]); z <- as.numeric(z[keep])
    n <- length(x)
    
    if (n < MIN_N) {
      skip_counter <- skip_counter + 1
      skipped[[skip_counter]] <- data.frame(
        cCRE_id = this_ccre, gene_id = gene, metabolite = metabolite,
        n = n, reason = paste0("n_below_", MIN_N))
      next
    }
    
    if (sd(x) == 0 || sd(y) == 0 || sd(z) == 0) {
      skip_counter <- skip_counter + 1
      skipped[[skip_counter]] <- data.frame(
        cCRE_id = this_ccre, gene_id = gene, metabolite = metabolite,
        n = n, reason = "zero_variance")
      next
    }
    
    #  Spearman partial correlation
    partial_corr <- tryCatch(
      ppcor::pcor.test(x, y, z, method = "spearman"),
      error = function(e) NULL
    )
    
    if (is.null(partial_corr) ||
        !is.finite(partial_corr$estimate) ||
        !is.finite(partial_corr$p.value)) {
      skip_counter <- skip_counter + 1
      skipped[[skip_counter]] <- data.frame(
        cCRE_id = this_ccre, gene_id = gene, metabolite = metabolite,
        n = n, reason = "invalid_partial_correlation")
      next
    }
    
    # numerical-stability diagnostic (rank-based determinant)
    rx <- rank(x); ry <- rank(y); rz <- rank(z)
    cor_matrix  <- cor(cbind(rx, ry, rz), method = "pearson")
    determinant <- det(cor_matrix)
    
    if (!is.finite(determinant)) {
      skip_counter <- skip_counter + 1
      skipped[[skip_counter]] <- data.frame(
        cCRE_id = this_ccre, gene_id = gene, metabolite = metabolite,
        n = n, reason = "nonfinite_determinant")
      next
    }
    
    unstable <- (abs(determinant) < 1e-8) ||
      (abs(as.numeric(partial_corr$estimate)) > 0.999)
    
    result_counter <- result_counter + 1
    results[[result_counter]] <- data.frame(
      cCRE_id              = this_ccre,
      gene_id              = gene,
      metabolite           = metabolite,
      n                    = n,
      spearman_partial_rho = as.numeric(partial_corr$estimate),
      partial_p            = as.numeric(partial_corr$p.value),
      determinant          = determinant,
      unstable             = unstable
    )
  }
  
  results <- bind_rows(results[seq_len(result_counter)])
  skipped <- bind_rows(skipped[seq_len(skip_counter)])
  
  write_csv(results, file.path(OUT_DIR, paste0("partial_correlation_", this_mark, ".csv")))
  write_csv(skipped, file.path(OUT_DIR, paste0("partial_correlation_", this_mark, "_skipped.csv")))
  
  cat(this_mark, " - valid: ", nrow(results),
      " | skipped: ", nrow(skipped),
      " | unstable(flagged): ", sum(results$unstable, na.rm = TRUE), "\n", sep = "")
}

# remove unstable rows,  BH-correct, both thresholds

generate_partial_fdr_outputs <- function(fdr_threshold) {
  
  suffix  <- ifelse(fdr_threshold == 0.05, "fdr05", "fdr10")
  sig_dir <- file.path(OUT_DIR, suffix)
  dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n==== Partial correlation FDR < ", fdr_threshold, " ====\n", sep = "")
  
  for (this_mark in top_marks) {
    
    in_file <- file.path(OUT_DIR, paste0("partial_correlation_", this_mark, ".csv"))
    results <- read_csv(in_file, show_col_types = FALSE)
    
    if (nrow(results) == 0) {
      cat(this_mark, " - no valid results, skipping\n", sep = "")
      next
    }
    
    # remove numerically unstable rows BEFORE correction
    clean_results <- results %>%
      dplyr::filter(!unstable, is.finite(partial_p)) %>%
      dplyr::mutate(partial_fdr = p.adjust(partial_p, method = "BH"))
    
    sig <- clean_results %>% dplyr::filter(partial_fdr < fdr_threshold)
    
    write_csv(clean_results,
              file.path(sig_dir, paste0("partial_correlation_", this_mark, "_all_clean_fdr.csv")))
    write_csv(sig,
              file.path(sig_dir, paste0("partial_correlation_", this_mark, "_significant.csv")))
    
    cat(this_mark,
        " - valid: ", nrow(results),
        " | after stability filter: ", nrow(clean_results),
        " | significant: ", nrow(sig), "\n", sep = "")
  }
}

generate_partial_fdr_outputs(0.05)
generate_partial_fdr_outputs(0.10)

#  summary table

summary_table <- data.frame(mark = top_marks)
summary_table[c("total_valid", "unstable", "valid_after_stability", "sig_fdr0.05", "sig_fdr0.1")]

for (i in seq_along(top_marks)) {
  
  this_mark <- top_marks[i]
  
  raw_file <- file.path(OUT_DIR, paste0("partial_correlation_", this_mark, ".csv"))
  raw <- read_csv(raw_file, show_col_types = FALSE)
  
  summary_table$total_valid[i]           <- nrow(raw)
  summary_table$unstable[i]              <- if (nrow(raw)) sum(raw$unstable, na.rm = TRUE) else 0
  summary_table$valid_after_stability[i] <- if (nrow(raw)) sum(!raw$unstable & is.finite(raw$partial_p)) else 0
  
  fdr05_file <- file.path(OUT_DIR, "fdr05", paste0("partial_correlation_", this_mark, "_significant.csv"))
  fdr10_file <- file.path(OUT_DIR, "fdr10", paste0("partial_correlation_", this_mark, "_significant.csv"))
  
  summary_table$sig_fdr05[i] <- if (file.exists(fdr05_file)) nrow(read_csv(fdr05_file, show_col_types = FALSE)) else NA
  summary_table$sig_fdr10[i] <- if (file.exists(fdr10_file)) nrow(read_csv(fdr10_file, show_col_types = FALSE)) else NA
}

print(summary_table)

write_csv(summary_table, file.path(OUT_DIR, "partial_correlation_summary_by_mark.csv"))
