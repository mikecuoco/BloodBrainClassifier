# load packages
library(reticulate) # interface with python
library(DESeq2) # differential expression

#' select DE genes from an RNA-seq count matrix
#' 
#' @param count.matrix a samples by genes matrix of non-negative numbers
#' @param conditions a set of sample labels for each row in count.matrix
#' @return DE genes with adjusted p-value < 0.05

select_DEgenes <- function(count.matrix, conditions){
  # for testing
  # start = import(module = "scripts.starter", as = "data")
  # count.matrix = start$X
  # conditions = start$y
  
  # create DESeq object
  coldata = data.frame(condition = as.factor(conditions), row.names = rownames(count.matrix))
  dds = DESeqDataSetFromMatrix(countData = round(t(count.matrix)), colData = coldata, design = ~ condition)

  # filter
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # run
  dds = DESeq(dds)
  res = results(dds)
  
  # return DE genes
  de.genes = na.omit(rownames(res)[res$padj < 0.05])
  return(as.vector(de.genes))
}