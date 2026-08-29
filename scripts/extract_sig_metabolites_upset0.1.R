combos <- read_csv(file.path(RESULTS_DIR, "FDR_0.1", "metabolites_per_combination_FDR0.1.csv"), show_col_types = FALSE)
combos %>%
      filter(n_marks == 5 |
       marks == "H3K4me3 + H3K27ac + H3K36me3 + H3K9me3" |
       marks == "H3K4me3 + H3K27me3 + H3K36me3 + H3K9me3") %>%
write_csv(file.path(RESULTS_DIR, "FDR_0.1", "three_arrow_metabolite_lists.csv"))

combos %>% filter(n_marks == 5) %>% pull(metabolite)
             
combos %>% filter(marks == "H3K4me3 + H3K27ac + H3K36me3 + H3K9me3") %>% pull(metabolite)

combos %>% filter(marks == "H3K4me3 + H3K27me3 + H3K36me3 + H3K9me3") %>% pull(metabolite)

#-------------
library(dplyr)
library(readr)

IN_DIR  <- file.path(RESULTS_DIR, "FDR_0.1")
OUT_DIR <- file.path(IN_DIR, "metabolite_lists")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

combos <- read_csv(file.path(IN_DIR, "metabolites_per_combination_FDR0.1.csv"),
                   show_col_types = FALSE)

# table 1: all five marks
table1 <- combos %>%
  filter(n_marks == 5) %>%
  dplyr::select(marks, metabolite)

cat("\n=== All 5 marks:", nrow(table1), "metabolites ===\n")
print(table1, n = Inf)

write_csv(table1, file.path(OUT_DIR, "metabolites_all5_marks.csv"))

# table 2: H3K4me3 + H3K27ac + H3K36me3 + H3K9me3
table2 <- combos %>%
  filter(marks == "H3K4me3 + H3K27ac + H3K36me3 + H3K9me3") %>%
  dplyr::select(marks, metabolite)

cat("\n=== H3K4me3 + H3K27ac + H3K36me3 + H3K9me3:", nrow(table2), "metabolites ===\n")
print(table2, n = Inf)

write_csv(table2, file.path(OUT_DIR, "metabolites_K4_K27ac_K36_K9.csv"))

# table 3: H3K4me3 + H3K27me3 + H3K36me3 + H3K9me3
table3 <- combos %>%
  filter(marks == "H3K4me3 + H3K27me3 + H3K36me3 + H3K9me3") %>%
  dplyr::select(marks, metabolite)

cat("\n=== H3K4me3 + H3K27me3 + H3K36me3 + H3K9me3:", nrow(table3), "metabolites ===\n")
print(table3, n = Inf)

write_csv(table3, file.path(OUT_DIR, "metabolites_K4_K27me3_K36_K9.csv"))

