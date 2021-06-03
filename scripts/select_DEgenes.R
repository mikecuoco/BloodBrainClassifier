#!/usr/bin/env Rscript

# load packages
library(reticulate) # interface with python
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2)) # differential expression
library(broom)

#' select DE genes from an RNA-seq count matrix
#' 
#' @param number of DE genes to keep (optional; defaults to 195)
#' @return top n DE genes with highest adjusted p-value

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  # defaullt values
  top_n = 100 
  method = "Ttest"
} else if (length(args)>0) {
  # default output file
  top_n = strtoi(args[1])
  method = args[2]
}

# load training data
start = import(module = "scripts.starter", as = "data")
count.matrix = start$X_train
conditions = start$y_train

if (method == "Ttest"){
  count.matrix = count.matrix %>% as.data.frame()
  # perform t-test on each gene
  test = map_df(count.matrix, function(.x){
    case = .x[conditions==1]
    control= .x[conditions==0]
    out = t.test(case,control) %>% tidy()
  })
  
  # wrangle
  res = test %>%
    mutate(gene = colnames(count.matrix)) %>%
    filter(!is.nan(p.value)) %>%
    mutate(padj = p.adjust(p.value, "fdr")) %>%
    column_to_rownames("gene")
  
  # return DE genes
  de.genes = as.vector(rownames(head(res[order(res$p.value),], top_n)))
  con = file("DEgenes.txt")
  writeLines(de.genes, con)
  close(con)
}

if (method ==  "DEseq"){
  # create DESeq object
  coldata = data.frame(condition = as.factor(conditions), row.names = rownames(count.matrix))
  dds = DESeqDataSetFromMatrix(countData = round(t(count.matrix)), colData = coldata, design = ~ condition)
  
  # remove genes with < 10 counts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # run
  dds = DESeq(dds)
  res = results(dds)
  
  # return DE genes
  #de.genes = as.vector(na.omit(rownames(res)[res$padj < 0.05])) # get genes with p-value less than 0.05, only
  de.genes = as.vector(rownames(head(res[order(res$padj),], top_n)))
  con = file("DEgenes.txt")
  writeLines(de.genes, con)
  close(con)
}



