library(jsonlite)
library(dplyr)
library(stringr)
library(readr)
metabolomic_ccle <- read_csv("/Users/emantaha/Downloads/CCLE_metabolomics_20190502.csv")

colnames(metabolomic_ccle)

encode_meta_data <- read.delim(
  "experiment_report_2026_3_12_9h_13m.tsv",
  stringsAsFactors = FALSE,
  skip = 1,  # skip first row
  header = TRUE,  
  fill = TRUE
)

encode_cell_lines <- encode_meta_data %>%
  filter(
    Assay.title == "Histone ChIP-seq",
    Biosample.classification == "cell line",
    Status == "released")

#===============================================================


# ENCODE cell lines that match with CCLE

meta_cell_lines <- metabolomic_ccle %>%
  distinct(CCLE_ID) %>%
  pull(CCLE_ID)


#remove special characters , convert strings to upper
clean_name <- function(x) {x %>%
    str_replace_all("[^a-zA-Z0-9]", "") %>%
    str_to_upper()
}

clean_ccle <- metabolomic_ccle %>%
  mutate(
    plain = str_remove(CCLE_ID, "_.*$"),
    clean = clean_name(plain)) %>%
  distinct(CCLE_ID, clean)

clean_encode <- encode_cell_lines %>%
  distinct(encode_name = Biosample.term.name) %>%
  mutate(clean = clean_name(encode_name))
  
matched_cell_lines <- inner_join(clean_encode, clean_ccle, by = "clean") %>%
  select(encode_name, ccle_id = CCLE_ID)

#missed Jurkat add manually 
matched_cell_lines <- matched_cell_lines %>%
  add_row(
    encode_name = "Jurkat, Clone E6-1", 
    ccle_id = "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"
  )

print(matched_cell_lines)

write_csv(matched_cell_lines, "results/encode_ccle_crossref.csv")

#=======================================================

# Frequent Histone marks & replicates


mark_rep_table <- encode_cell_lines %>%
  filter(Biosample.term.name %in% matched_cell_lines$encode_name) %>%
  group_by(cell_line = Biosample.term.name,
           mark = Target.of.assay) %>%
  summarise(n_experiments = n(),
            n_bio_reps = n_distinct(Biological.replicate),
            accessions = paste(Accession, collapse = " | ")
            ) %>%
  add_count(mark, name = "n_cell_lines_per_mark")%>%
  arrange(desc(n_cell_lines_per_mark), mark, cell_line)

# total marks per cell line

total_marks <- mark_rep_table %>%
  group_by(cell_line) %>%
  summarise(
    n_marks      = n_distinct(mark),
    marks        = paste(sort(unique(mark)), collapse = ", "),
    total_experiments = sum(n_experiments),
    avg_bio_reps = round(mean(n_bio_reps), 1)
  ) %>%
  arrange(desc(n_marks))
write_csv(mark_rep_table, "results/Histone_marks_frequency.csv")
write_csv(total_marks, "results/total_marks.csv")

most_freq_marks <- mark_rep_table %>%
  distinct(mark, n_cell_lines_per_mark) %>%
  arrange(desc(n_cell_lines_per_mark))

write_csv(most_freq_marks, "results/most_frequant_marks.csv")

#=========================================================





