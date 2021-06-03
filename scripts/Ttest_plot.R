# #!/usr/bin/env Rscript

# load packages
library(reticulate) # interface with python
suppressPackageStartupMessages(library(tidyverse))
library(broom)
# load training data
start = import(module = "scripts.starter", as = "data")
count.matrix = start$X_train
count.matrix = count.matrix[,colSums(count.matrix) > 0]
conditions = start$y_train
genes = start$genes %>% rownames_to_column("ensembl")

# perform t-test on each gene
count.matrix = count.matrix %>% as.data.frame()
test = map_df(count.matrix, function(.x){
  case = .x[conditions==1]
  control= .x[conditions==0]
  out = t.test(case,control) %>% tidy()
})

# wrangle
res = test %>%
  mutate(ensembl = colnames(count.matrix)) %>%
  filter(!is.nan(p.value)) %>%
  mutate(padj = p.adjust(p.value, "fdr")) %>%
  left_join(genes, by = "ensembl") 

theme_set(theme_minimal())
rbind(top_n(res,30, statistic), top_n(res,-10, statistic)) %>%
  ggplot(aes(y = reorder(gene_name,statistic), x = statistic, fill = statistic < 0)) +
  geom_bar(stat="identity") +
  labs(y = NULL, x = "t-statistic") +
  theme(legend.position = "none")

ggsave("Tstatistic_bars.png", height = 8, width = 5)

theme_set(theme_classic())
res %>%
  ggplot(aes(x=statistic)) +
  geom_histogram(fill = "gray", color = "black") +
  labs(x = "t-statistic")

ggsave("Tstatistic_dist.png", height = 4, width = 6)
