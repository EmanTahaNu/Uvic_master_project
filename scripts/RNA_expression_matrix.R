library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(httr)
library(jsonlite)
library(purrr)

# folders
BASE_DIR     <- "/Users/emantaha/Desktop/thesis"
RNA_DIR      <- file.path(BASE_DIR, "rna_quant")
MATRICES_DIR <- file.path(BASE_DIR, "matrices")
dir.create(RNA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MATRICES_DIR, showWarnings = FALSE, recursive = TRUE)

#=====================================================
# load ENCODE report
#=====================================================

encode_RNA_data <- read.delim(file.path(BASE_DIR, "last_RNA_encode_data.tsv"),
                              skip = 1, header = TRUE, sep = "\t",
                              quote = "", fill = TRUE, stringsAsFactors = FALSE)

#=====================================================
# filter: total RNA-seq, Barbara + Thomas lab only
#=====================================================

rna_sheet <- encode_RNA_data %>%
  filter(Assay.title == "total RNA-seq",
         Lab %in% c("Barbara Wold, Caltech", "Thomas Gingeras, CSHL")) %>%
  mutate(clean_name = str_replace_all(Biosample.term.name, "[^a-zA-Z0-9]", "")) %>%
  select(Biosample.term.name, Accession, Assay.title, Lab,
         Biological.replicate, clean_name, Files)

write.csv(rna_sheet, file.path(BASE_DIR, "results/RNA_sheet.csv"), row.names = FALSE)

#=====================================================
# keep only my 27 cell lines
#=====================================================

matched_names <- matched_cell_lines$encode_name

rna_sheet <- rna_sheet %>%
  filter(Biosample.term.name %in% matched_names)

cat("missing cell lines:\n")
print(setdiff(matched_names, unique(rna_sheet$Biosample.term.name)))

write.csv(rna_sheet, file.path(BASE_DIR, "results/RNA_sheet_filtered.csv"), row.names = FALSE)

#=====================================================
# expand Files column: one row per file accession
#=====================================================

rna_files <- rna_sheet %>%
  mutate(file_id = strsplit(as.character(Files), ",")) %>%
  unnest(file_id) %>%
  mutate(file_id = trimws(file_id),
         file_accession = str_extract(file_id, "ENCFF[A-Z0-9]+")) %>%
  filter(!is.na(file_accession))

#=====================================================
# ask ENCODE for real download link + file type
#=====================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

get_file_info <- function(acc) {
  r <- httr::GET(paste0("https://www.encodeproject.org/files/", acc, "/?format=json"),
                 httr::add_headers(Accept = "application/json"))
  
  if (httr::status_code(r) != 200) {
    return(tibble(file_accession = acc, status = NA, output_type = NA, href = NA))
  }
  
  x <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"))
  
  tibble(file_accession = acc,
         status      = x$status %||% NA,
         output_type = x$output_type %||% NA,
         href        = x$href %||% NA)
}

file_info <- purrr::map_dfr(unique(rna_files$file_accession), get_file_info)

rna_files <- rna_files %>%
  left_join(file_info, by = "file_accession") %>%
  filter(!is.na(status), !is.na(output_type),
         status == "released", output_type == "gene quantifications") %>%
  mutate(file_url   = paste0("https://www.encodeproject.org", href),
         quant_file = file.path(RNA_DIR, paste0(clean_name, "_RNAseq_", Accession, "_",
                                                file_accession, ".tsv")))

cat("files to download:", nrow(rna_files), "\n")

#=====================================================
# download
#=====================================================

for (i in seq_len(nrow(rna_files))) {
  url <- rna_files$file_url[i]
  f   <- rna_files$quant_file[i]
  
  if (is.na(url)) {
    cat("no url:", rna_files$file_accession[i], "\n")
    next
  }
  
  download.file(url, f, mode = "wb", quiet = FALSE)
  cat("downloaded:", rna_files$file_accession[i], "\n")
  Sys.sleep(0.5)
}

#=====================================================
# read files, check gene_id + TPM exist
#=====================================================

exp_list <- list()

for (i in seq_len(nrow(rna_files))) {
  f <- rna_files$quant_file[i]
  if (!file.exists(f)) next
  
  first_line <- readLines(f, n = 1, warn = FALSE)
  if (length(first_line) == 0 || grepl("<!DOCTYPE html>", first_line, ignore.case = TRUE)) {
    cat("bad file, skipping:", basename(f), "\n")
    next
  }
  
  dat <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(dat) || !all(c("gene_id", "TPM") %in% names(dat))) {
    cat("missing gene_id/TPM:", basename(f), "\n")
    next
  }
  
  exp_list[[length(exp_list) + 1]] <- dat %>%
    select(gene_id, TPM) %>%
    mutate(cell_line = rna_files$Biosample.term.name[i])
}

all_exp <- bind_rows(exp_list)
cat("cell lines loaded:\n")
print(unique(all_exp$cell_line))

#=====================================================
# strip gene version so same gene matches across cell lines
#=====================================================

all_exp$gene_id <- str_remove(all_exp$gene_id, "\\..*")

#=====================================================
# collapse replicates using median
#=====================================================

collapsed <- all_exp %>%
  group_by(gene_id, cell_line) %>%
  summarise(TPM = median(TPM, na.rm = TRUE), .groups = "drop")

#=====================================================
# build matrix, keep ENSG genes only
#=====================================================

mat <- collapsed %>%
  pivot_wider(names_from = cell_line, values_from = TPM) %>%
  filter(grepl("^ENSG", gene_id))

cell_cols <- setdiff(colnames(mat), "gene_id")
cat("NA count:", sum(is.na(mat[, cell_cols])), "\n")

#=====================================================
# log2(TPM+1) version
#=====================================================

mat_log <- mat
mat_log[, cell_cols] <- log2(as.matrix(mat_log[, cell_cols]) + 1)

#=====================================================
# round both
#=====================================================

mat[, cell_cols]     <- round(mat[, cell_cols], 4)
mat_log[, cell_cols] <- round(mat_log[, cell_cols], 4)

#=====================================================
# save
#=====================================================

saveRDS(mat,     file.path(MATRICES_DIR, "expression_matrix2_TPM.rds"))
saveRDS(mat_log, file.path(MATRICES_DIR, "expression_matrix2_log2TPM.rds"))
write_tsv(mat,     file.path(MATRICES_DIR, "expression_matrix2_TPM.tsv"))
write_tsv(mat_log, file.path(MATRICES_DIR, "expression_matrix2_log2TPM.tsv"))

