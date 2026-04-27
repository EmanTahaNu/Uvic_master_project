library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)
library(httr)
library(jsonlite)

#==================================================
dir.create("/Users/emantaha/Desktop/thesis/rna_quant",
           showWarnings = FALSE, recursive = TRUE)
BASE_DIR     <- "/Users/emantaha/Desktop/thesis"
RNA_DIR      <- file.path(BASE_DIR, "rna_quant")
MATRICES_DIR <- file.path(BASE_DIR, "matrices")
dir.create(RNA_DIR,      showWarnings = FALSE)
dir.create(MATRICES_DIR, showWarnings = FALSE)


encode_RNA_data <- read.delim("experiment_report_2026_4_22_5h_0m.tsv"
  ,
  stringsAsFactors = FALSE,
  skip = 1,  # skip first row
  header = TRUE,  
  fill = TRUE
)

#finding all RNA-seq experiments for cell lines

rna_sheet <- encode_RNA_data %>%
  filter( Assay.title  %in% c("total RNA-seq", "polyA plus RNA-seq"),
          Biosample.classification == "cell line",
          Status == "released",
          Biosample.term.name %in% matched_cell_lines$encode_name) %>%
  mutate(
    clean_name = str_replace_all(Biosample.term.name, "[^a-zA-Z0-9]", ""),
    quant_file = file.path(RNA_DIR,
                           paste0(clean_name, "_RNAseq_", Accession,".tsv"))
  ) %>% 
  dplyr::select(cell_line = Biosample.term.name,
         assay = Assay.title,
         accession = Accession,
         clean_name, quant_file)
print(rna_sheet %>% count(cell_line, name = "n_exp") %>% count(n_exp))

write.csv(rna_sheet, "results/RNA_sheet.csv")

#=================================================
# Download Tsv files from ENCODE(V29) with Ensembl IDs

get_tsv_url <- function(accession) {
  url <- paste0("https://www.encodeproject.org/experiments/",
                accession, "/@@embedded?format=json")
  resp <- tryCatch(GET(url, timeout(30)), error = function(e) NULL)
  if(is.null(resp) || http_error(resp)) return(NA_character_)
  
  files <- fromJSON(content(resp, "text", encoding = "UTF-8"),
                    simplifyVector = TRUE)$files
  
  if(is.null(files) || nrow(files) == 0) return(NA_character_)
  
  match <- files %>%
    filter(file_format == "tsv",
           output_type == "gene quantifications",
           assembly == "GRCh38",
           status == "released")
  
  v29 <- match %>% filter(grepl("V29|v29", genome_annotation))
  
  if(nrow(v29) > 0) {
    return(paste0("https://www.encodeproject.org", v29$href[1]))
  }
  
  
  if (nrow(match) > 0) {
    return(paste0("https://www.encodeproject.org", match$href[1]))
 
  }
  return(NA_character_)
  
}
  
#--------------------------------------------

# Download loop

for (i in seq_len(nrow(rna_sheet))) {
  
  f <- rna_sheet$quant_file[i]
  label <- paste(rna_sheet$cell_line[i], rna_sheet$accession[i])
  
  if(file.exists(f)) {cat ("skip (exists):", rna_sheet$cell_line[i], "\n") ;
    next}
  
  url <- get_tsv_url(rna_sheet$accession[i])
  
if(is.na(url)) { cat("no url found", label, "\n"); next }
  
  system(paste0('curl -L -o "', f, '" "', url, '"'))
  Sys.sleep(1)
}

#------------------------------------------  
# Read files 

rna_sheet <- rna_sheet %>% filter(file.exists(quant_file))

all_exp_list <- vector("list", nrow(rna_sheet))

for (i in seq_len(nrow(rna_sheet))) {
  all_exp_list[[i]] <- read_tsv(rna_sheet$quant_file[i], show_col_types = FALSE) %>%
    dplyr::select(gene_id, TPM) %>%
    mutate(cell_line = rna_sheet$cell_line[i],
           accession = rna_sheet$accession[i])}  

all_exp <- bind_rows(all_exp_list)

print(unique(all_exp$cell_line))

#----------------------------------------------------
#collapse replicates using median
collapsed <- all_exp %>%
  group_by(gene_id, cell_line) %>%
  summarise(TPM = median(TPM, na.rm = TRUE), .groups = "drop")

#reshape matrix

mat <- collapsed %>%
  pivot_wider(names_from = cell_line, values_from = TPM)

cell_cols <- setdiff(colnames(mat), "gene_id")

head(cell_cols)

# Replace NA with 0

mat <- mat %>%
  mutate(across(all_of(cell_cols), ~ replace_na(.x, 0)))



na_count <- sum(is.na(mat[, cell_cols]))

# log2 + 1 transform

mat_log   <- mat
mat_log[, cell_cols] <- log2(as.matrix(mat_log[, cell_cols]) + 1)

mat <- mat %>% filter(grepl("^ENSG", gene_id))
mat_log <- mat_log %>% filter(grepl("ENSG", gene_id))
head(mat$gene_id)

# round TPM to 4 decimals
mat[, cell_cols]  <- round(mat[, cell_cols], 4)

mat_log[, cell_cols]  <- round(mat_log[, cell_cols], 4)

#save matrix
saveRDS(mat,     file.path(MATRICES_DIR, "expression_matrix2_TPM.rds"))
saveRDS(mat_log, file.path(MATRICES_DIR, "expression_matrix2_log2TPM.rds"))
write_tsv(mat,   file.path(MATRICES_DIR, "expression_matrix2_TPM.tsv"))
write_tsv(mat,   file.path(MATRICES_DIR, "expression_matrix2_log2TPM.tsv"))







  
  
  
  
  
  
  
  
  
  
  
  










