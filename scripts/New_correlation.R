
source("/Volumes/eman/thesis/thesis_paths.R")
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(cowplot)

top_marks <- c("H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3")


CORR_RERUN_DIR <- file.path(RESULTS_DIR, "pls_correlation_full_rerun")

dir.create(CORR_RERUN_DIR, recursive = TRUE)

# ================================================================
#load inputs

ccre_bed <- read_tsv(
  file.path(CCRE_DIR, "GRCh38-cCREs.bed"),
  col_names = c("chr", "start", "end", "cCRE_id", "score", "category"),
  show_col_types = FALSE
)

pls_ids <- ccre_bed %>% dplyr::filter(category == "PLS") %>% dplyr::pull(cCRE_id)

met <- read_csv(file.path(MATRICES_DIR, "metabolite_matrix_aligned.csv"),
                show_col_types = FALSE) %>%
  tibble::column_to_rownames("metabolite")

# ================================================================
#run correlation for every mark

for (this_mark in top_marks) {
  
  out_file <- file.path(CORR_RERUN_DIR, paste0("PLS_correlation_", this_mark, ".csv"))
  
  cat("Computing correlations for", this_mark, "\n")
  
  h_mat <- readRDS(file.path(ALL_MATRICES_DIR,
                             paste0("histone_matrix_", this_mark, "_aligned.rds")))
  
  h_pls <- h_mat %>% dplyr::filter(cCRE_id %in% pls_ids)
  cat(this_mark, ": ", nrow(h_pls), "PLS cCREs\n", sep = "")
  
  common_cells <- intersect(setdiff(colnames(h_pls), c("cCRE_id", "category")),
                            colnames(met))
  cat(this_mark, ": ", length(common_cells), "common cell lines\n", sep = "")
  
  h_vals <- as.matrix(h_pls[, common_cells])
  m_vals <- as.matrix(met[, common_cells])
  
  results_list <- vector("list", nrow(h_vals))
  
  for (i in seq_len(nrow(h_vals))) {
    if (i %% 1000 == 0) cat("  Progress:", i, "/", nrow(h_vals), "\n")
    
    histone_vec <- as.numeric(h_vals[i, ])
    
    row_results <- data.frame(
      cCRE_id    = h_pls$cCRE_id[i],
      metabolite = rownames(m_vals),
      pearson_r  = NA_real_,
      pearson_p  = NA_real_,
      spearman_r = NA_real_,
      spearman_p = NA_real_
    )
    
    for (j in seq_len(nrow(m_vals))) {
      met_vec <- as.numeric(m_vals[j, ])
      
      if (var(histone_vec, na.rm = TRUE) == 0 ||
          var(met_vec, na.rm = TRUE) == 0) next
      
      p_test <- cor.test(histone_vec, met_vec, method = "pearson", use = "pairwise.complete.obs")
      s_test <- cor.test(histone_vec, met_vec, method = "spearman", use = "pairwise.complete.obs")
      
      row_results$pearson_r[j]  <- p_test$estimate
      row_results$pearson_p[j]  <- p_test$p.value
      row_results$spearman_r[j] <- s_test$estimate
      row_results$spearman_p[j] <- s_test$p.value
    }
    
    results_list[[i]] <- row_results
  }
  
  results <- dplyr::bind_rows(results_list) %>%
    dplyr::mutate(
      mark         = this_mark,
      pearson_fdr  = p.adjust(pearson_p, method = "BH"),
      spearman_fdr = p.adjust(spearman_p, method = "BH")
    )
  
  write_csv(results, out_file)
}

# ================================================================
#threshold outputs, run once per FDR level

generate_fdr_outputs <- function(fdr_threshold) {
  
  suffix  <- ifelse(fdr_threshold == 0.05, "fdr05", "fdr0.1")
  out_dir <- file.path(CORR_RERUN_DIR, paste0("correlation_", suffix))
  
  dir.create(out_dir, recursive = TRUE)
  
  all_mark_results <- list()
  
  # per-mark significant counts, density plots 
  for (this_mark in top_marks) {
    
    results <- read_csv(file.path(CORR_RERUN_DIR, paste0("PLS_correlation_", this_mark, ".csv")),
                        show_col_types = FALSE)
    all_mark_results[[this_mark]] <- results
    
    n_sig_pearson  <- sum(results$pearson_fdr  < fdr_threshold, na.rm = TRUE)
    n_sig_spearman <- sum(results$spearman_fdr < fdr_threshold, na.rm = TRUE)
    cat(this_mark, "- sig Pearson:", n_sig_pearson, " sig Spearman:", n_sig_spearman, "\n")
    
    # density: all pairs
    all_pairs <- results %>%
      dplyr::filter(!is.na(pearson_r), !is.na(spearman_r)) %>%
      tidyr::pivot_longer(cols = c(pearson_r, spearman_r), names_to = "method", values_to = "r") %>%
      dplyr::mutate(method = ifelse(method == "pearson_r", "Pearson", "Spearman"))
    
    ggplot(all_pairs, aes(x = r, fill = method, color = method)) +
      geom_density(alpha = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      scale_fill_manual(values  = c("Pearson" = "steelblue", "Spearman" = "tomato")) +
      scale_color_manual(values = c("Pearson" = "steelblue", "Spearman" = "tomato")) +
      labs(title = paste("All pairs —", this_mark, "at PLS cCREs"),
           subtitle = paste("n =", nrow(results), "pairs | FDR <", fdr_threshold),
           x = "Correlation coefficient (r)", y = "Density") +
      theme_classic() + theme(legend.position = "top")
    ggsave(file.path(out_dir, paste0("density_allpairs_", this_mark, ".png")), width = 7, height = 5, dpi = 300)
    
    # density: significant pairs only
    sig_pairs <- results %>%
      dplyr::filter(spearman_fdr < fdr_threshold | pearson_fdr < fdr_threshold) %>%
      dplyr::filter(!is.na(pearson_r), !is.na(spearman_r)) %>%
      tidyr::pivot_longer(cols = c(pearson_r, spearman_r), names_to = "method", values_to = "r") %>%
      dplyr::mutate(method = ifelse(method == "pearson_r", "Pearson", "Spearman"))
    
    if (nrow(sig_pairs) == 0) {
      cat("  No significant pairs at FDR <", fdr_threshold, "for", this_mark, "\n")
    } else {
      ggplot(sig_pairs, aes(x = r, fill = method, color = method)) +
        geom_density(alpha = 0.4) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        scale_fill_manual(values  = c("Pearson" = "steelblue", "Spearman" = "tomato")) +
        scale_color_manual(values = c("Pearson" = "steelblue", "Spearman" = "tomato")) +
        labs(title = paste("Significant pairs only —", this_mark, "at PLS cCREs"),
             subtitle = paste("n =", nrow(sig_pairs) / 2, "significant pairs | FDR <", fdr_threshold),
             x = "Correlation coefficient (r)", y = "Density") +
        theme_classic() + theme(legend.position = "top")
      ggsave(file.path(out_dir, paste0("density_sigpairs_", this_mark, ".png")), width = 7, height = 5, dpi = 300)
    }
  }
  
 
  
 
 
  # combined density figure 
  make_plot <- function(df, title) {
    df %>%
      dplyr::filter(!is.na(pearson_r), !is.na(spearman_r)) %>%
      tidyr::pivot_longer(all_of(c("pearson_r", "spearman_r")), names_to = "method", values_to = "r") %>%
      dplyr::mutate(method = ifelse(method == "pearson_r", "Pearson", "Spearman")) %>%
      ggplot(aes(x = r, fill = method, color = method)) +
      geom_density(alpha = 0.4, linewidth = 0.6) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      scale_fill_manual(values  = c("Pearson" = "#E64B35", "Spearman" = "#4DBBD5")) +
      scale_color_manual(values = c("Pearson" = "#E64B35", "Spearman" = "#4DBBD5")) +
      scale_x_continuous(limits = c(-1, 1), breaks = c(-0.5, 0, 0.5)) +
      labs(title = title, x = "r", y = "Density") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
            legend.position = "none", axis.title = element_text(size = 8), axis.text = element_text(size = 7))
  }
  
  legend <- cowplot::get_legend(
    make_plot(all_mark_results[[1]], "") + theme(legend.position = "bottom", legend.title = element_blank())
  )
  
  row_all <- cowplot::plot_grid(
    plotlist = lapply(top_marks, function(m) make_plot(all_mark_results[[m]], m)), nrow = 2)
  
  row_sig <- cowplot::plot_grid(
    plotlist = lapply(top_marks, function(m) {
      df <- all_mark_results[[m]] %>% dplyr::filter(spearman_fdr < fdr_threshold)
      if (nrow(df) == 0) return(ggplot() + theme_void() +
                                  annotate("text", x = 0.5, y = 0.5, label = "No sig. pairs", color = "gray50", size = 3))
      make_plot(df, m)
    }), nrow = 2)
  
  final <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label("A  All pairs", fontface = "bold", size = 11, x = 0.02, hjust = 0),
    row_all,
    cowplot::ggdraw() + cowplot::draw_label(paste0("B  Significant pairs only (FDR < ", fdr_threshold, ")"),
                                            fontface = "bold", size = 11, x = 0.02, hjust = 0),
    row_sig,
    legend,
    ncol = 1, rel_heights = c(0.06, 1, 0.06, 1, 0.1)
  )
  
  ggsave(file.path(out_dir, "combined_density_publication.png"), final, width = 14, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, "combined_density_publication.pdf"), final, width = 14, height = 7, bg = "white")
  
  # unique significant pairs / cCREs / metabolites per mark
  for (this_mark in top_marks) {
    results <- all_mark_results[[this_mark]]
    sig <- results %>% dplyr::filter(spearman_fdr < fdr_threshold | pearson_fdr < fdr_threshold)
    cat(this_mark, "- total sig pairs:", nrow(sig),
        " unique cCREs:", dplyr::n_distinct(sig$cCRE_id),
        " unique metabolites:", dplyr::n_distinct(sig$metabolite), "\n")
  }
  
  #  significant cCREs per mark bar plot 
  mark_counts <- data.frame(
    mark = top_marks,
    n_ccres = sapply(top_marks, function(m) {
      all_mark_results[[m]] %>%
        dplyr::filter(spearman_fdr < fdr_threshold | pearson_fdr < fdr_threshold) %>%
        dplyr::pull(cCRE_id) %>% unique() %>% length()
    })
  ) %>% dplyr::mutate(mark = factor(mark, levels = top_marks))
  
  ggplot(mark_counts, aes(x = mark, y = n_ccres, fill = mark)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = scales::comma(n_ccres)), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("H3K4me3" = "#E64B35", "H3K27ac" = "#4DBBD5", "H3K27me3" = "#7B2D8B",
                                 "H3K36me3" = "#F0A500", "H3K9me3" = "#00897B")) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Unique significant cCREs per histone mark",
         subtitle = paste("PLS promoters — Spearman FDR <", fdr_threshold),
         x = "Histone mark", y = "Number of unique significant cCREs") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "gray40", hjust = 0.5),
          legend.position = "none")
  
  ggsave(file.path(out_dir, "cCRE_per_mark_barplot.png"), width = 7, height = 5, dpi = 300, bg = "white")
  
}

# ================================================================
# run for both thresholds

generate_fdr_outputs(0.05)
generate_fdr_outputs(0.1)

