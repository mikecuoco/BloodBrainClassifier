#!/usr/bin/env Rscript

# load packages
library(reticulate) # interface with python
suppressPackageStartupMessages(library(DESeq2)) # differential expression

#' select DE genes from an RNA-seq count matrix
#' 
#' @param number of DE genes to keep (optional; defaults to 195)
#' @return top n DE genes with highest adjusted p-value

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  top_n = 195 # genes with p-values less than 0.05
} else if (length(args)==1) {
  # default output file
  top_n = strtoi(args[1])
}

# load training data
start = import(module = "scripts.starter", as = "data")
count.matrix = start$X_train
conditions = start$y_train

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
de.genes = as.vector(na.omit(rownames(res)[res$padj < 0.05])) # get genes with p-value less than 0.05, only
de.genes = as.vector(rownames(head(res[order(res$padj),], top_n)))
cat(de.genes,sep=",")

