library(dplyr)
library(readr)
library(UpSetR)

top_marks <- c("H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3")
FDR_CUTOFF <- 0.05
FULL_DIR <- file.path(RESULTS_DIR, "pls_correlation_full_rerun")
# significant metabolite names per mark
met_sets <- lapply(setNames(top_marks, top_marks), function(m) {
  results <- read_csv(file.path(FULL_DIR,
                                paste0("PLS_correlation_", m, ".csv")),
                                      show_col_types = FALSE)
  results %>%
    filter(spearman_fdr < FDR_CUTOFF | pearson_fdr < FDR_CUTOFF) %>%
    pull(metabolite) %>%
    unique()
})

sapply(met_sets, length)

# binary matrix for upset
all_mets <- unique(unlist(met_sets))
metabolite_matrix <- as.data.frame(
  sapply(met_sets, function(s) as.integer(all_mets %in% s))
)
rownames(metabolite_matrix) <- all_mets

png(file.path(IN_DIR, "upset_metabolites_FDR0.05.png"), width = 1200, height = 800)
upset(metabolite_matrix,
      sets            = top_marks,
      order.by        = "freq",
      sets.bar.color  = "red",
      main.bar.color  = "red",
      text.scale      = 1.3,
      mainbar.y.label = "Number of Metabolites",
      sets.x.label    = "Metabolites per Mark",
      mb.ratio        = c(0.6, 0.4))
dev.off()

# exclusive metabolite combinations (same thresholds/sets as the plot above)
all_combos <- data.frame()

for (n in 5:2) {
  for (combo in combn(top_marks, n, simplify = FALSE)) {
    common <- Reduce(intersect, met_sets[combo])
    for (o in setdiff(top_marks, combo)) common <- setdiff(common, met_sets[[o]])
    if (length(common) > 0)
      all_combos <- bind_rows(all_combos,
                              data.frame(n_marks    = n,
                                         marks      = paste(combo, collapse = " + "),
                                         metabolite = common))
  }
}

print(all_combos %>% arrange(desc(n_marks)))

write_csv(all_combos, file.path(RESULTS_DIR, "metabolites_per_combination_FDR0.05.csv"))


#------------------------------------------


FDR_CUTOFF <- 0.1
IN_DIR <- file.path(RESULTS_DIR, "FDR_0.1")
FULL_DIR <- file.path(RESULTS_DIR, "pls_correlation_full_rerun")
# significant metabolite names per mark
met_sets <- lapply(setNames(top_marks, top_marks), function(m) {
  results <- read_csv(file.path(FULL_DIR,
                                paste0("PLS_correlation_", m, ".csv")),
                      show_col_types = FALSE)
  results %>%
    filter(spearman_fdr < FDR_CUTOFF | pearson_fdr < FDR_CUTOFF) %>%
    pull(metabolite) %>%
    unique()
})

sapply(met_sets, length)

# binary matrix for upset
all_mets <- unique(unlist(met_sets))
metabolite_matrix <- as.data.frame(
  sapply(met_sets, function(s) as.integer(all_mets %in% s))
)
rownames(metabolite_matrix) <- all_mets

png(file.path(IN_DIR, "upset_metabolites_FDR0.1.png"), width = 1200, height = 800)
upset(metabolite_matrix,
      sets            = top_marks,
      order.by        = "freq",
      sets.bar.color  = "red",
      main.bar.color  = "red",
      text.scale      = 1.3,
      mainbar.y.label = "Number of Metabolites",
      sets.x.label    = "Metabolites per Mark",
      mb.ratio        = c(0.6, 0.4))
dev.off()

# exclusive metabolite combinations (same thresholds/sets as the plot above)
all_combos <- data.frame()

for (n in 5:2) {
  for (combo in combn(top_marks, n, simplify = FALSE)) {
    common <- Reduce(intersect, met_sets[combo])
    for (o in setdiff(top_marks, combo)) common <- setdiff(common, met_sets[[o]])
    if (length(common) > 0)
      all_combos <- bind_rows(all_combos,
                              data.frame(n_marks    = n,
                                         marks      = paste(combo, collapse = " + "),
                                         metabolite = common))
  }
}

print(all_combos %>% arrange(desc(n_marks)))

write_csv(all_combos, file.path(IN_DIR, "metabolites_per_combination_FDR0.1.csv"))




