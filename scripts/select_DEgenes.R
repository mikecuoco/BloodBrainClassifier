#!/usr/bin/env Rscript

# load packages
library(reticulate) # interface with python
suppressPackageStartupMessages(library(DESeq2)) # differential expression

#' select DE genes from an RNA-seq count matrix
#' 
#' @param count.matrix a samples by genes matrix of non-negative numbers
#' @param conditions a vector of sample labels for each row in count.matrix
#' @return DE genes with adjusted p-value < 0.05

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
de.genes = as.vector(na.omit(rownames(res)[res$padj < 0.05]))
# start$data = start$data[,c(de.genes,"group")]
cat(de.genes,sep=",")

