
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

MACDB   <- "/Volumes/eman/thesis/downloads.metabolite.txt"
SIG_DIR <- "/Volumes/eman/thesis/results/FDR_0.1/significant"
OUT_DIR <- "/Volumes/eman/thesis/results/enrichment_shared_metabolites"
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

marks <- c("H3K4me3","H3K27ac","H3K27me3","H3K36me3","H3K9me3")

# rebuild the 5-mark metabolite list directly 
all_sig <- lapply(marks, function(m){
  f <- file.path(SIG_DIR, paste0("PLS_correlation_", m, "_significant.csv"))
  if (file.exists(f)) read_csv(f, show_col_types=FALSE) %>% mutate(mark=m) else NULL
}) %>% bind_rows()

met_marks <- all_sig %>%
  distinct(metabolite, mark) %>%
  group_by(metabolite) %>%
  summarise(n_marks = n_distinct(mark),
            which_marks = paste(sort(mark), collapse="+"),
            .groups="drop")

mets_5 <- met_marks %>% filter(n_marks == 5) %>%
  dplyr::select(metabolite, which_marks)

cat("Pan-mark metabolites:", nrow(mets_5), "\n")
print(mets_5$metabolite)

# save it this time
write_csv(mets_5, file.path(OUT_DIR, "metabolites_5mark.csv"))



#  load MACdb
mac <- read_tsv(MACDB, show_col_types=FALSE)
mac$pubchem_CID <- suppressWarnings(as.character(mac$pubchem_CID))

#  normalised name match
norm <- function(x) x %>% tolower() %>% str_replace_all("[^a-z0-9]","")
mac_names <- unique(mac$original_metabolite_name)
mac_norm  <- norm(mac_names)

match_tbl <- tibble(your_metabolite = mets_5$metabolite) %>%
  mutate(macdb_name = mac_names[match(norm(my_metabolite), mac_norm)],
         matched    = !is.na(macdb_name))

cat("\nMatched to MACdb:", sum(match_tbl$matched),
    "of", nrow(match_tbl), "\n")
cat("\nNOT matched:\n")
print(match_tbl %>% filter(!matched) %>% pull(my_metabolite))

#  pull MACdb records for matched metabolites 
key  <- match_tbl %>% filter(matched) %>% dplyr::select(macdb_name, my_metabolite)
hits <- mac %>%
  filter(original_metabolite_name %in% key$macdb_name) %>%
  left_join(key, by=c("original_metabolite_name"="macdb_name"))

#  summary 
summary_5 <- hits %>%
  mutate(p = suppressWarnings(as.numeric(`case_control_p-value`))) %>%
  group_by(my_metabolite) %>%
  summarise(
    n_records   = n(),
    n_cohorts   = n_distinct(Cohort_id),
    n_sig_p05   = sum(p < 0.05, na.rm=TRUE),
    prop_sig    = round(n_sig_p05 / n_records, 2),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_p05))

print(summary_5, n=30)

write_csv(summary_5, file.path(OUT_DIR, "macdb_5mark_metabolites.csv"))

#plot: all 24 metabolites
# merge unmatched as n=0 so all 24 appear
all_24 <- mets_5 %>%
  left_join(summary_5, by=c("metabolite"="your_metabolite")) %>%
  mutate(n_sig_p05 = ifelse(is.na(n_sig_p05), 0, n_sig_p05),
         n_records = ifelse(is.na(n_records), 0, n_records),
         in_macdb  = n_records > 0)

p <- ggplot(all_24, aes(x = reorder(metabolite, n_sig_p05),
                        y = n_sig_p05,
                        fill = in_macdb)) +
  geom_col() +
  geom_text(aes(label = ifelse(n_sig_p05 > 0, n_sig_p05, "not in MACdb")),
            hjust = -0.1, size = 3.2) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE"  = "#2c7fb8",
                               "FALSE" = "#d3d3d3"),
                    labels = c("TRUE"  = "In MACdb",
                               "FALSE" = "Not found"),
                    name = NULL) +
  labs(title = "Cancer evidence for the 24 pan-mark metabolites (MACdb)",
       subtitle = "Number of studies with significant case/control difference (p<0.05)",
       x = NULL, y = "Significant cancer studies") +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face="bold"),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "macdb_5mark_check.png"), p, width=9, height=8, dpi=150)
cat("\nSaved: macdb_5mark_check.png\n")

print(summary_5 %>%
        filter(my_metabolite %in% c("2-hydroxyglutarate","alpha-ketoglutarate",
                                      "citrate","isocitrate","oxalate",
                                      "phosphocreatine","pyroglutamic acid",
                                      "5-adenosylhomocysteine","thiamine")))




match_4 <- tibble(your_metabolite = mets_4) %>%
  mutate(macdb_name = mac_names[match(norm(your_metabolite), mac_norm)],
         matched    = !is.na(macdb_name))

cat("\n====== MACdb match — 4-mark metabolites ======\n")
cat("Matched:", sum(match_4$matched), "of", nrow(match_4), "\n")
cat("\nNot matched:\n")
print(match_4 %>% filter(!matched) %>% pull(your_metabolite))

# pull records 
key4  <- match_4 %>% filter(matched) %>% dplyr::select(macdb_name, your_metabolite)
hits4 <- mac %>%
  filter(original_metabolite_name %in% key4$macdb_name) %>%
  left_join(key4, by=c("original_metabolite_name"="macdb_name"))

# summary 
summary_4 <- hits4 %>%
  mutate(p = suppressWarnings(as.numeric(`case_control_p-value`))) %>%
  group_by(your_metabolite) %>%
  summarise(
    n_records  = n(),
    n_cohorts  = n_distinct(Cohort_id),
    n_sig_p05  = sum(p < 0.05, na.rm=TRUE),
    prop_sig   = round(n_sig_p05 / n_records, 2),
    .groups    = "drop"
  ) %>%
  arrange(desc(n_sig_p05))

print(summary_4, n=50)

write_csv(summary_4, file.path(OUT_DIR, "macdb_4mark_metabolites.csv"))

#plot 
all_4 <- tibble(metabolite = mets_4) %>%
  left_join(summary_4, by=c("metabolite"="your_metabolite")) %>%
  mutate(n_sig_p05 = ifelse(is.na(n_sig_p05), 0, n_sig_p05),
         n_records  = ifelse(is.na(n_records),  0, n_records),
         in_macdb   = n_records > 0)

p_mac4 <- ggplot(all_4, aes(x = reorder(metabolite, n_sig_p05),
                            y = n_sig_p05,
                            fill = in_macdb)) +
  geom_col() +
  geom_text(aes(label = ifelse(n_sig_p05 > 0,
                               as.character(n_sig_p05),
                               "not in MACdb")),
            hjust = -0.1, size = 3.1) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE"  = "#bd0026",
                               "FALSE" = "#d3d3d3"),
                    labels = c("TRUE"  = "In MACdb",
                               "FALSE" = "Not found"),
                    name   = NULL) +
  labs(title    = "Cancer evidence for 4-mark shared metabolites (MACdb)",
       subtitle = "Number of studies with significant case/control difference (p<0.05)",
       x = NULL, y = "Significant cancer studies") +
  theme_minimal(base_size = 12) +
  theme(plot.title     = element_text(face="bold"),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "macdb_4mark_check.png"),
       p_mac4, width=9, height=8, dpi=150)

sum5 <- read_csv(file.path(OUT_DIR, "macdb_5mark_metabolites.csv"),
                 show_col_types=FALSE) %>%
  mutate(group = "5 marks (pan-mark)")

sum4 <- summary_4 %>% mutate(group = "4 marks")

combined <- bind_rows(sum5, sum4) %>%
  filter(n_sig_p05 > 0) %>%
  arrange(desc(n_sig_p05))

print(combined %>% dplyr::select(your_metabolite, group, n_records,
                                 n_cohorts, n_sig_p05, prop_sig), n=50)

write_csv(combined, file.path(OUT_DIR, "macdb_combined_4_5mark.csv"))
