
library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

SIG_DIR <- "/Volumes/eman/thesis/results/FDR_0.1/significant"
ASSIGN  <- "/Volumes/eman/thesis/results/gene_assignment_v29/PLS_gene_assignment_v29_pc.txt"
OUT_DIR <- "/Volumes/eman/thesis/results/enrichment_shared_metabolites"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

marks <- c("H3K4me3","H3K27ac","H3K27me3","H3K36me3","H3K9me3")

# load all significant pairs 
all_sig <- lapply(marks, function(m){
  f <- file.path(SIG_DIR, paste0("PLS_correlation_", m, "_significant.csv"))
  if (file.exists(f)) read_csv(f, show_col_types=FALSE) %>% mutate(mark=m) else NULL
}) %>% bind_rows()

# ---- 2. how many marks does each metabolite reach? ----
met_marks <- all_sig %>%
  distinct(metabolite, mark) %>%
  group_by(metabolite) %>%
  summarise(n_marks      = n_distinct(mark),
            which_marks  = paste(sort(mark), collapse="+"),
            .groups = "drop")

# the two sets of interest 
mets_5mark <- met_marks %>% filter(n_marks == 5) %>% pull(metabolite)
mets_4mark <- met_marks %>% filter(n_marks == 4) %>% pull(metabolite)

cat("Metabolites shared across ALL 5 marks:", length(mets_5mark), "\n")
cat(paste(" ", sort(mets_5mark)), sep="\n")
cat("\nMetabolites shared across exactly 4 marks:", length(mets_4mark), "\n")
cat(paste(" ", sort(mets_4mark)), sep="\n")

#  load gene assignment & background universe 
assign <- read_tsv(ASSIGN, show_col_types=FALSE)
assign$gene_id <- sub("\\..*$", "", assign$gene_id)

universe_symbols <- unique(assign$gene_name)
universe_entrez  <- unique(bitr(universe_symbols, "SYMBOL", "ENTREZID",
                                OrgDb = org.Hs.eg.db)$ENTREZID)
cat("\nBackground universe:", length(universe_entrez), "genes\n")

# enrichment function 
enrich_metabolite_set <- function(metabolite_set, label, colour="#2c7fb8") {
  

  # cCREs associated with these metabolites (across any mark)
  ccres <- all_sig %>%
    filter(metabolite %in% metabolite_set) %>%
    pull(cCRE_id) %>% unique()
  cat("cCREs:", length(ccres), "\n")
  
  # genes at those cCREs
  genes <- assign %>%
    filter(cCRE_id %in% ccres) %>%
    pull(gene_name) %>% unique()
  cat("Genes:", length(genes), "\n")
  
  if (length(genes) < 10) {
    cat("Too few genes for enrichment — skipping\n")
    return(NULL)
  }
  
  # symbol -> Entrez
  fg <- unique(bitr(genes, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
  cat("Foreground Entrez IDs:", length(fg), "\n")
  
  # GO BP
  ego <- enrichGO(gene          = fg,
                  universe      = universe_entrez,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.10,
                  readable      = TRUE)
  
  n_terms <- if (is.null(ego)) 0 else nrow(as.data.frame(ego))
  cat("GO BP terms (FDR<0.05/q<0.10):", n_terms, "\n")
  
  if (n_terms > 0) {
    # save table
    write_csv(as.data.frame(ego),
              file.path(OUT_DIR, paste0("GO_BP_", label, ".csv")))
    
    # dotplot — top 20 terms
    p <- dotplot(ego, showCategory = 20) +
      scale_colour_gradient(low = "grey80", high = colour) +
      ggtitle(paste("GO Biological Process —", gsub("_"," ", label)),
              subtitle = paste(length(metabolite_set), "shared metabolites |",
                               length(ccres), "cCREs |",
                               length(genes), "genes")) +
      theme_minimal(base_size = 12) +
      theme(plot.title    = element_text(face="bold", size=13),
            plot.subtitle = element_text(size=10, colour="grey40"),
            axis.text.y   = element_text(size=9))
    ggsave(file.path(OUT_DIR, paste0("GO_BP_", label, ".png")),
           p, width=9, height=8, dpi=150)
    cat("Saved:", paste0("GO_BP_", label, ".png"), "\n")
    
    # print top 10 terms
    cat("\nTop 10 terms:\n")
    print(as.data.frame(ego) %>%
            select(Description, GeneRatio, p.adjust) %>%
            head(10))
  }
  
  return(ego)
}

#  run enrichment
ego_5 <- enrich_metabolite_set(mets_5mark, "5mark_shared", colour = "#08519c")
ego_4 <- enrich_metabolite_set(mets_4mark, "4mark_shared", colour = "#bd0026")

write_csv(
  met_marks %>% filter(n_marks == 5) %>% dplyr::select(metabolite, which_marks),
  file.path(OUT_DIR, "metabolites_5mark.csv"))
write_csv(
  met_marks %>% filter(n_marks == 4) %>% dplyr::select(metabolite, which_marks),
  file.path(OUT_DIR, "metabolites_4mark.csv"))

cat("\n\n= SUMMARY =\n")
cat("5-mark metabolites:", length(mets_5mark), "\n")
cat("4-mark metabolites:", length(mets_4mark), "\n")



library(dplyr); library(readr); library(clusterProfiler)
library(org.Hs.eg.db); library(ggplot2)

SIG_DIR <- "/Volumes/eman/thesis/results/FDR_0.1/significant"
ASSIGN  <- "/Volumes/eman/thesis/results/gene_assignment_v29/PLS_gene_assignment_v29_pc.txt"
OUT_DIR <- "/Volumes/eman/thesis/results/enrichment_shared_metabolites"
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

marks <- c("H3K4me3","H3K27ac","H3K27me3","H3K36me3","H3K9me3")

# rebuild significant pairs
all_sig <- lapply(marks, function(m){
  f <- file.path(SIG_DIR, paste0("PLS_correlation_", m, "_significant.csv"))
  if (file.exists(f)) read_csv(f, show_col_types=FALSE) %>% mutate(mark=m) else NULL
}) %>% bind_rows()

met_marks <- all_sig %>%
  distinct(metabolite, mark) %>%
  group_by(metabolite) %>%
  summarise(n_marks     = n_distinct(mark),
            which_marks = paste(sort(mark), collapse="+"),
            .groups     = "drop")

# the 4-mark set
mets_4 <- met_marks %>% filter(n_marks == 4) %>% pull(metabolite)
cat("4-mark metabolites:", length(mets_4), "\n")
print(sort(mets_4))

# gene assignment + background
assign <- read_tsv(ASSIGN, show_col_types=FALSE)
assign$gene_id <- sub("\\..*$", "", assign$gene_id)
universe_entrez <- unique(bitr(unique(assign$gene_name),
                               "SYMBOL", "ENTREZID",
                               OrgDb=org.Hs.eg.db)$ENTREZID)
cat("Background:", length(universe_entrez), "genes\n")

# cCREs and genes for the 4-mark metabolites
ccres_4 <- all_sig %>% filter(metabolite %in% mets_4) %>% pull(cCRE_id) %>% unique()
genes_4  <- assign %>% filter(cCRE_id %in% ccres_4) %>% pull(gene_name) %>% unique()
cat("cCREs:", length(ccres_4), "| Genes:", length(genes_4), "\n")

fg_4 <- unique(bitr(genes_4, "SYMBOL", "ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID)

cat("Foreground Entrez IDs:", length(fg_4), "\n")

# GO enrichment
ego_4 <- enrichGO(gene          = fg_4,
                  universe      = universe_entrez,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.10,
                  readable      = TRUE)

n4 <- if (is.null(ego_4)) 0 else nrow(as.data.frame(ego_4))
cat("GO BP terms (4-mark set):", n4, "\n")

if (n4 > 0) {
  # save table
  write_csv(as.data.frame(ego_4),
            file.path(OUT_DIR, "GO_BP_4mark_shared.csv"))
  
  # dotplot
  p <- dotplot(ego_4, showCategory=20) +
    scale_colour_gradient(low="grey80", high="#bd0026") +
    ggtitle("GO Biological Process — 4-mark shared metabolites",
            subtitle = paste(length(mets_4), "metabolites |",
                             length(ccres_4), "cCREs |",
                             length(genes_4), "genes")) +
    theme_minimal(base_size=12) +
    theme(plot.title=element_text(face="bold", size=13),
          plot.subtitle=element_text(size=10, colour="grey40"))
  ggsave(file.path(OUT_DIR, "GO_BP_4mark_shared.png"),
         p, width=9, height=8, dpi=150)
  cat("Saved: GO_BP_4mark_shared.png\n")
  
  # print top terms safely
  cat("\nTop 15 GO terms:\n")
  res <- as.data.frame(ego_4)
  print(res[1:min(15,nrow(res)), c("Description","GeneRatio","p.adjust")])
}

# save metabolite list
write_csv(
  met_marks %>% filter(n_marks==4) %>% dplyr::select(metabolite, which_marks),
  file.path(OUT_DIR, "metabolites_4mark.csv"))

cat("\nDone. Outputs in:", OUT_DIR, "\n")
