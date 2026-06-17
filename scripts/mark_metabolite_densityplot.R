library(dplyr)
library(ggplot2)
library(readr)
library(cowplot)
library(purrr)

RESULTS_DIR <- "/Users/emantaha/Desktop/thesis/results"
marks <- c("H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3")
geodist_dir <- file.path(RESULTS_DIR, "geodist_direction")
density_dir <- file.path(RESULTS_DIR, "density_single")
combined_dir <- file.path(RESULTS_DIR, "combined")

dir.create(single_dir, showWarnings = FALSE)

dir.create(combined_dir, showWarnings = FALSE)


for (mark in marks) {
  for (direction in c("pos", "neg")) {
    
    files <- list.files(geodist_dir,
                        pattern = paste0("^dist_", dir, "_", mark, "_.*\\.txt$"),
                        full.names = TRUE)
    
   if (length(files) == 0) next
    
    color <- "steelblue"
    if (direction == "neg") color <- "red"
    
    plots <- list()
    for (f in files) {
      dat <- read_tsv(f, show_col_types = FALSE)
      
      if (nrow(dat) < 2) next
      if (length(unique(dat$distance_bp)) < 2) next
      
      dat$distance_kb <- dat$distance_bp / 1000
      
      name <- tools::file_path_sans_ext(basename(f))
      title <- sub(paste0("^dist_", direction, "_", mark, "_"), "", name)
      
      p <- ggplot(dat, aes(distance_kb)) +
        geom_density(fill = color, alpha = 0.05) +
        scale_x_log10() +
        labs(title = title, x = "Distance (kb)", y = "Density") +
        theme_classic()
      
      ggsave(file.path(single_dir, paste0(name, ".png")), p, width = 6, height = 4)
      
      plots[[length(plots) + 1]] <- p
      
    }
      
    if (length(plots) == 0 ) next
    n_rows <- ceiling(length(plots) / 4)
    
    title_label <- "Positive"
    if (direction == "neg") title_label <- "Negative"
    
    banner <- ggdraw() + draw_label(paste(mark, "-", title_label))
    grid <- plot_grid(plotlist = plots, ncol = 4)
    final <- plot_grid(banner, grid, ncol = 1)
    
    out_file <- file.path(combined_dir, paste0("combined_", direction, "_", mark, ".png"))
    ggsave(out_file, final)
    
  }
  
}
    
    
    
    
    
    
    