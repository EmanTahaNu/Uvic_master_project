library(dplyr)
library(readr)


# 1. Metabolite x cell line matrix

metabolite_matrix <- metabolomic_ccle %>%
  filter(CCLE_ID %in% matched_cell_lines$ccle_id) %>%
  select(-DepMap_ID) %>%
  tibble::column_to_rownames("CCLE_ID") %>%
  t() %>%
  as.data.frame()

write_csv(metabolite_matrix %>% tibble::rownames_to_column("metabolite"), "matrices/metabolite_matrix.csv")
 
#===========================================

dir.create("bigwigs", showWarnings = FALSE )
dir.create("ccre_signal", showWarnings = FALSE)
dir.create("matrices", showWarnings = FALSE)
dir.create("cCRE", showWarnings = FALSE)

# Top Marks:

top_marks <- c("H3K4me3", "H3K27ac", "H3K27me3","H3K36me3", "H3K9me3" )
ccre_bed <- "cCRE/GRCh38-cCREs.bed"

# Construct sample_sheet

sample_sheet <- mark_rep_table %>%
  filter(mark %in% top_marks) %>%
  mutate(
    cell_line_clean = str_replace_all(cell_line, "[^a-zA-Z0-9]", ""),
    bigwig_file = paste0("bigwigs/", cell_line_clean, "_", mark, ".bigwig"),
    signal_file = paste0("ccre_signal/", cell_line_clean, "_", mark, "_signal.tab")
  )


write_csv(sample_sheet, "results/sample_sheet_csv")


sample_sheet <- sample_sheet %>%
  mutate(
    cell_line_clean = str_replace_all(cell_line, "[^a-zA-Z0-9]", ""),
    bigwig_file  = paste0("/Users/emantaha/Desktop/thesis/bigwigs/",
                          cell_line_clean, "_", mark, ".bigWig"),
    signal_file  = paste0("/Users/emantaha/Desktop/thesis/ccre_signal/",
                          cell_line_clean, "_", mark, "_signal.tab")
  )
#=========================================================

#Get bigwig URL from ENCODE API

get_bigwig_url <- function(accession) {
  url <- paste0("https://www.encodeproject.org/experiments/", accession, "/@@embedded?format=json")

  response <- tryCatch(
    GET(url, add_headers(Accept = "application/json"), timeout(30)),
    error = function(e) NULL
  )
  
  if(is.null(response) || http_error(response)) return(NA_character_)
  
  files <- fromJSON(content(response, as = "text", encoding = "UTF-8"),
                    simplifyVector = FALSE)[["files"]]

  if (is.null(files) || length(files) == 0) return(NA_character_)
  
  for (f in files) {
    if (isTRUE(f$file_format == "bigWig")                   &&
        isTRUE(f$output_type == "fold change over control") &&
        isTRUE(f$assembly    == "GRCh38")                   &&
        isTRUE(f$status      == "released")) {
  
  return(paste0("https://www.encodeproject.org", f$href))
  
   }
  }
  return(NA_character_)
}


#===================================================
# Download BigWig files 

for (i in seq_len(nrow(sample_sheet))) {
  bw_file <- sample_sheet$bigwig_file[i]
  accession <- sample_sheet$accession[i]
  label <- paste(sample_sheet$cell_line[i], sample_sheet$mark[i])
  
  if (file.exists(bw_file)) { cat("SKIP", label, "\n"); next }
  
  bw_url <- get_bigwig_url(accession)
  if(is.na(bw_url))  {cat("NOT FOUND", label, "\n") ;next}
  
  cat("Downloading", i, "/", nrow(sample_sheet), ":", label, "\n")
  
  system(paste0('curl -L --progress-bar -o "', bw_file, '" "',  bw_url, '"'))
  Sys.sleep(1)
  
}



#===========================================

# Signal computation

ccre_bed <- "/Users/emantaha/Desktop/thesis/cCRE/GRCh38-cCREs-3col.bed"

for (i in seq_len(nrow(sample_sheet))) {
  
  bw_file <- sample_sheet$bigwig_file[i]
  signal_file <- sample_sheet$signal_file[i]
  label <- paste(sample_sheet$cell_line[i], sample_sheet$mark[i])
  
  if (!file.exists(bw_file))  {cat("SKIP no bigwig:", label, "\n"); next}
  if (file.exists(signal_file)) {cat("SKIP done:", label, "\n"); next}
  
  system(paste("bigWigAverageOverBed", bw_file, ccre_bed, signal_file ))
  

}

#Emans-MacBook-Air:~ emantaha$ chmod +x bigWigAverageOverBed
#Emans-MacBook-Air:~ emantaha$ sudo mv bigWigAverageOverBed /usr/local/bin/
#  Password:
 
#Emans-MacBook-Air:~ emantaha$ cd Desktop/thesis/cCRE/
#  Emans-MacBook-Air:cCRE emantaha$ ls
#GRCh38-cCREs.bed
#Emans-MacBook-Air:cCRE emantaha$ head -5 GRCh38-cCREs.bed 
#chr1	10033	10250	EH38D4327497	EH38E2776516	pELS
#chr1	10385	10713	EH38D4327498	EH38E2776517	pELS
#chr1	16097	16381	EH38D6144701	EH38E3951272	CA-CTCF
#chr1	17343	17642	EH38D6144702	EH38E3951273	CA-TF
#chr1	29320	29517	EH38D6144703	EH38E3951274	CA


#===================================================

# cCRE x cell line matrices

read_signal <- function(file_path, cell_line_name) {
  read_tsv(file_path, 
           col_names = c("cCRE_id", "size", "covered", "sum",
                         "mean0", "mean"),
           col_types = cols(cCRE_id = col_character(),
                            .default = col_double()),
           show_col_types = FALSE) %>%
    select(cCRE_id, !!cell_line_name := mean)
}

top_marks <- c("H3K4me3", "H3K27ac", "H3K27me3","H3K36me3", "H3K9me3" )


for (mark in top_marks) {
  
  rows <- sample_sheet %>%
    filter(mark == !!mark, file.exists(signal_file))
  
  if (nrow(rows) == 0) {cat( "no signal files found", mark,"\n"); next }
  
  mat <- read_signal(rows$signal_file[1], rows$cell_line[1])
  
  for (j in 2:nrow(rows)) {
    mat <- full_join(mat, read_signal(rows$signal_file[j], rows$cell_line[j]),
                     by = "cCRE_id")
    
  }
  
  #remove rows with zero variance
  
  columns <- setdiff(colnames(mat), "cCRE_id")
  rm_v <- apply(mat[ ,columns], 1, var, na.rm =TRUE)
  mat <- mat[!is.na(rm_v) & rm_v > 0, ]
  
saveRDS(mat, paste0("/Users/emantaha/Desktop/thesis/matrices/histone_matrix_", mark, ".rds"))
write_tsv(mat, gzfile(paste0("/Users/emantaha/Desktop/thesis/matrices/histone_matrix_", mark, ".tsv.gz")))
          
  
}













