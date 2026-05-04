library(ggplot2)

# read the cCRE bed file
ccre_bed <- read_tsv("/Users/emantaha/Desktop/thesis/cCRE/GRCh38-cCREs.bed",
                     col_names = c("chr", "start", "end", "cCRE_id", "score", "category"),
                     show_col_types = FALSE)

head(ccre_bed)

print(table(ccre_bed$category))


# count cCREs per category
category_counts <- ccre_bed %>%
  count(category) %>%
  arrange(desc(n))

# plot
ggplot(category_counts, aes(x = reorder(category, -n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = scales::comma(n)), vjust = -0.5, size = 3.5) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Number of cCREs per Category",
    x     = "Category",
    y     = "Number of cCREs"
  ) +
  theme_classic() 

# save
ggsave(file.path(RESULTS_DIR, "cCRE_category_barplot.png"),
       width = 7, height = 5, dpi = 300)








