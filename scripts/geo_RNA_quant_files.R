library(GEOquery)
library(dplyr)

RESULTS_DIR <- "/Users/emantaha/Desktop/thesis/results"
INPUT_FILE <- "/Users/emantaha/Desktop/thesis/untreated_geo_deduplicateded.tsv"

clean_data <- read.delim(INPUT_FILE, stringsAsFactors = FALSE)
gse_list <- unique(clean_data$gse_accession)
cat("total unique GSE accessions to check:", length(gse_list), "\n")

results <- data.frame()

for (i in seq_along(gse_list)) {
  
  gse <- gse_list[i]
  cat(i, "/", length(gse_list), " checking", gse, "\n")
  
  file_list <- tryCatch(
    getGEOSuppFiles(GEO = gse, fetch_files = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(file_list) || nrow(file_list) == 0) {
    results <- rbind(results, data.frame(
      gse_accession = gse, supp_file_names = NA, status = "no_supp_file"
    ))
    next
  }
  
  fnames <- file_list$fname
  quant_like <- fnames[grepl("count|tpm|fpkm|rsem|quant|expr|matrix", fnames, ignore.case = TRUE)]
  
  results <- rbind(results, data.frame(
    gse_accession = gse,
    supp_file_names = paste(fnames, collapse = " | "),
    status = ifelse(length(quant_like) > 0, "has_quant_supp_file", "has_supp_file_unclear_type")
  ))
}

results <- results %>%
  left_join(
    clean_data %>% group_by(gse_accession) %>% summarise(cell_lines = paste(sort(unique(cell_line)), collapse = ";")),
    by = "gse_accession"
  ) %>%
  relocate(cell_lines, .after = gse_accession)

write.csv(results, file.path(RESULTS_DIR, "geo_supp_file_check2.csv"), row.names = FALSE)

print(table(results$status))

