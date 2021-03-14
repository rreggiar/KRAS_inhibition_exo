#!/user/bin/env Rscript

## rreggiar@ucsc.edu

## testing job functionality in Rstudio for running big chunks in background

## load in libs -- these are apparently not included in the environment
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggsci)
  library(ggpubr)
  library(ggsignif)
  library(ggthemes) 
  library(ggrepel)
  library(ggforce)
  library(ggpmisc)
  library(ggprism)
  library(patchwork)
  library(extrafont)
  library(ggrepel)
  library(pheatmap)
  library(viridis)
  library(patchwork)
  library(fgsea)
  library(caret)
  library(msigdbr)
  library(matrixStats)
  library(NMF)
  library(rsample)
  library(vip)
  library(ROCR)
  library(glmnet)
  library(mlbench)
  library(tximeta)
  library(SummarizedExperiment)
  library(rjson)
  library(HDF5Array)
  library(DESeq2)
  library(ggpmisc)
  library(sva)
})
suppressMessages(loadfonts())

## read in toil rsem counts (log21p)
toil.expected.counts <- read.delim(file.path(output.data.dir, 'tcga.data', 'TcgaTargetGtex_gene_expected_count.gz')) %>%
  select(sample, contains('TCGA')) 
## convert from wide to long
toil.expected.counts %>%
  rename(sample = 'ensg') %>% 
  gather(sample, count, -ensg) -> toil.expected.counts.tidy
## subset phenotype df to only the samples of interest
luad.phenotypes.merge %>% filter(KRAS %in% c('p.G12C', 'WT_normal')) -> luad.kras.merge
## merge in phenotype columns
toil.expected.counts.tidy %>% 
  mutate(sample = gsub(sample, pattern = '[.]', replacement = '-')) %>% 
  merge(luad.kras.merge, by = 'sample') %>% 
  distinct() -> toil.expected.counts.tidy.merge
## merge in gene names for interpretability
toil.expected.counts.final <- 
  toil.expected.counts.tidy.merge %>% 
  merge(gen.24.names, by = 'ensg') %>%
  select(-ensg)
## generate DESeq compatible colData
toil.expected.counts.final %>% 
  select(sample, sample_type, KRAS) %>% 
  distinct() -> toil.expected.counts.coldata
## convert back to wide, undo log transform 
### RERWIP across(where()) seems awfully slow..
toil.expected.counts.final %>% 
  select(sample, gene, count) %>% 
  pivot_wider(names_from = gene, values_from = count, values_fn = mean) %>% 
  mutate(across(where(is.numeric), ~round((2^.) - 1))) -> toil.expected.counts.unlog.matrix
## transpose matrix to match DESeq expectations
toil.expected.counts.unlog.matrix %>% column_to_rownames('sample') %>% 
  t() %>% as.data.frame() -> toil.expected.counts.matrix.t
## create DESeq object with colData and counts
toil.dds <- DESeqDataSetFromMatrix(countData = toil.expected.counts.matrix.t,
                                   colData = toil.expected.counts.coldata,
                                   design = ~KRAS)
## filter by suggested abundance threshold
toil.dds.keep <- rowSums(counts(toil.dds)) >= 10
toil.dds <- toil.dds[toil.dds.keep,]
## DESeq
toil.dds <- DESeq(toil.dds)
## available comparisons
print(resultsNames(toil.dds))
## chosen comparison
print(tail(resultsNames(toil.dds), n=1))
## MA plot to assess health of model
plotMA(lfcShrink(toil.dds, coef = tail(resultsNames(toil.dds), n=1)))
## shrink, parse, and export DE results
kras_wt.vs.g12c_de <- lfcShrink(toil.dds, coef = tail(resultsNames(toil.dds), n=1)) %>% 
  as.data.frame() %>% 
  filter(padj <= 0.05) %>% 
  rownames_to_column('gene') %>% 
  select(gene, baseMean, log2FoldChange, padj) %>% 
  mutate(log2FoldChange = -log2FoldChange)