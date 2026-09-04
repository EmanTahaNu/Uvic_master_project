# ================================================================
# GO + KEGG enrichment on V29 gene assignment , both FDR thresholds
# Background: all PLS cCREs assigned to a protein-coding gene (V29)


BASE   <- "/Volumes/eman/thesis/results"
ASSIGN <- file.path(BASE, "gene_assignment_v29", "PLS_gene_assignment_v29_pc.txt")
OUT    <- file.path(BASE, "enrichment_v29")
dir.create(OUT, showWarnings = FALSE)

marks <- c("H3K4me3","H3K27ac","H3K27me3","H3K36me3","H3K9me3")

# two runs: label -> (input dir, threshold)
runs <- list(
  fdr05 = list(dir = file.path(BASE, "significant_correlation"), cut = 0.05),
  fdr0.1 = list(dir = file.path(BASE, "FDR_0.1", "significant"),   cut = 0.1)
)

assign <- read_tsv(ASSIGN, show_col_types = FALSE)
universe_symbols <- unique(assign$gene_name)
universe_entrez  <- unique(bitr(universe_symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)
cat("Background universe genes (Entrez):", length(universe_entrez), "\n\n")

for (lbl in names(runs)) {
  SIG_DIR <- runs[[lbl]]$dir
  FDR_CUT <- runs[[lbl]]$cut
  
  
  for (mk in marks) {
    f <- file.path(SIG_DIR, paste0("PLS_correlation_", mk, "_significant.csv"))
    if (!file.exists(f)) { cat(mk, "- no file, skip\n"); next }
    
    sig <- read_csv(f, show_col_types = FALSE)
    if (nrow(sig) == 0) { cat(mk, "- 0 pairs, skip\n\n"); next }
    if ("spearman_fdr" %in% names(sig)) sig <- sig %>% filter(spearman_fdr < FDR_CUT)
    if (nrow(sig) == 0) { cat(mk, "- 0 after filter, skip\n\n"); next }
    
    sig_genes <- assign %>%
      filter(cCRE_id %in% unique(sig$cCRE_id)) %>%
      pull(gene_name) %>% unique()
    
    cat(mk, "- sig cCREs:", length(unique(sig$cCRE_id)),
        "| genes:", length(sig_genes), "\n")
    if (length(sig_genes) < 5) { cat("  too few genes, skip\n\n"); next }
    
    fg <- unique(bitr(sig_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)
    
    ego <- enrichGO(gene = fg, universe = universe_entrez, OrgDb = org.Hs.eg.db,
                    ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE)
    ekegg <- enrichKEGG(gene = fg, universe = universe_entrez, organism = "hsa",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05)
    
    nGO <- if (is.null(ego))  0 else nrow(as.data.frame(ego))
    nKE <- if (is.null(ekegg)) 0 else nrow(as.data.frame(ekegg))
    cat("  GO BP:", nGO, "| KEGG:", nKE, "\n")
    
    if (nGO > 0) {
      write_csv(as.data.frame(ego), file.path(OUT, paste0("GO_BP_", mk, "_", lbl, ".csv")))
      ggsave(file.path(OUT, paste0("GO_BP_", mk, "_", lbl, ".png")),
             dotplot(ego, showCategory = 15) + ggtitle(paste("GO BP —", mk, lbl)),
             width = 8, height = 6, dpi = 150)
    }
    if (nKE > 0) {
      write_csv(as.data.frame(ekegg), file.path(OUT, paste0("KEGG_", mk, "_", lbl, ".csv")))
      ggsave(file.path(OUT, paste0("KEGG_", mk, "_", lbl, ".png")),
             dotplot(ekegg, showCategory = 15) + ggtitle(paste("KEGG —", mk, lbl)),
             width = 8, height = 6, dpi = 150)
    }
    cat("\n")
  }
}
