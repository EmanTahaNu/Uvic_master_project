library(dplyr)
library(readr)
library(ggplot2)

top_marks <- c("H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3")

FULL_DIR <- file.path(RESULTS_DIR, "pls_correlation_full_rerun")
OUT_DIR  <- file.path(RESULTS_DIR, "distance_distribution")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

MAX_DIST   <- 250e6   # 250 million bp
FDR_CUTOFF <- 0.1

# ---- load cCRE coordinates ----
ccre_bed <- read_tsv(
  "/Volumes/eman/thesis/cCRE/GRCh38-cCREs.bed",
  col_names = c("chr", "start", "end", "cCRE_id", "score", "category"),
  show_col_types = FALSE
) %>%
  mutate(midpoint = (start + end) / 2)

# ---- collect significant cCREs across all marks ----
sig_ccres <- character()
for (this_mark in top_marks) {
  res <- read_csv(file.path(FULL_DIR, paste0("PLS_correlation_", this_mark, ".csv")),
                  show_col_types = FALSE)
  sig <- res %>%
    filter(spearman_fdr < FDR_CUTOFF | pearson_fdr < FDR_CUTOFF) %>%
    pull(cCRE_id) %>%
    unique()
  sig_ccres <- union(sig_ccres, sig)
}

cat("Total unique significant cCREs (all marks):", length(sig_ccres), "\n")

# ---- keep only their coordinates ----
sig_coords <- ccre_bed %>% filter(cCRE_id %in% sig_ccres)
cat("Matched to coordinates:", nrow(sig_coords), "\n")

# ---- pairwise distances within each chromosome ----
all_dist <- data.frame()

for (this_chr in unique(sig_coords$chr)) {
  
  chr_ccres <- sig_coords %>% filter(chr == this_chr) %>% arrange(midpoint)
  n <- nrow(chr_ccres)
  if (n < 2) next
  
  # all unique pairs on this chromosome
  pairs <- combn(n, 2)
  d <- abs(chr_ccres$midpoint[pairs[1, ]] - chr_ccres$midpoint[pairs[2, ]])
  
  all_dist <- bind_rows(all_dist, data.frame(chr = this_chr, distance = d))
}

cat("Total same-chromosome pairs:", nrow(all_dist), "\n")

# ---- cap at 250 Mbp ----
all_dist <- all_dist %>% filter(distance <= MAX_DIST)
cat("Pairs within 250 Mbp:", nrow(all_dist), "\n")

write_csv(all_dist, file.path(OUT_DIR, "ccre_pairwise_distances.csv"))

# ---- distribution plot ----
ggplot(all_dist, aes(x = distance / 1e6)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +
  labs(title = "Pairwise distance between significant cCREs (same chromosome)",
       subtitle = paste0(nrow(all_dist), " pairs, capped at 250 Mbp"),
       x = "Distance (Mbp)", y = "Number of pairs") +
  theme_classic()

ggsave(file.path(OUT_DIR, "ccre_distance_distribution.png"),
       width = 8, height = 5, dpi = 300, bg = "white")

# ---- zoomed view for short-range structure (<= 5 Mbp) ----
short <- all_dist %>% filter(distance <= 5e6)

ggplot(short, aes(x = distance / 1e6)) +
  geom_histogram(bins = 60, fill = "tomato", color = "white") +
  labs(title = "Short-range distances (<= 5 Mbp)",
       subtitle = paste0(nrow(short), " pairs"),
       x = "Distance (Mbp)", y = "Number of pairs") +
  theme_classic()

ggsave(file.path(OUT_DIR, "ccre_distance_distribution_short.png"),
       width = 8, height = 5, dpi = 300, bg = "white")