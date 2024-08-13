#####
# author: "Hongfei Xu" (hongfeifly)
# title: "**Data pretreatment for Go enrichment analysis**"
# output:
#   pdf_document: default
#   html_document: default
#   fig_width: 4
#   fig_height: 4
#   fontsize: 20pt
# descriptions: transform log2Foldchange between castes for each tissue to 0-1.
#####

## setup work directory
setwd("rna_seq_analysis/Go")

library(tidyverse)

###########################################################################################
## Import gene expression data of each tissue
b_gene <- read.csv("data/brain_gene_between_castes.csv")
colnames(b_gene)[1] <- "gene_id"
b_gene <- b_gene[-c(1:4),]

f_gene <- read.csv("data/fatbody_gene_between_castes.csv")
colnames(f_gene)[1] <- "gene_id"
f_gene <- f_gene[-c(1:4),]

o_gene <- read.csv("data/ovary_gene_between_castes.csv")
colnames(o_gene)[1] <- "gene_id"
o_gene <- o_gene[-c(1:4),]

s_gene <- read.csv("data/spermatheca_gene_between_castes.csv")
colnames(s_gene)[1] <- "gene_id"
s_gene <- s_gene[-c(1:4),]

## define two functions to normalize log2Foldchange to 0-1/1-0
min_max_normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

max_min_normalize <- function(x) {
  return((x - max(x)) / (min(x) - max(x)))
}

## normalize log2FC value for each data
b_queen_score <- b_gene %>%
  mutate(allScore = min_max_normalize(b_gene$log2FoldChange)) %>%
  dplyr::select(gene_id, allScore, log2FoldChange)

b_worker_score <- b_gene %>%
  mutate(allScore = max_min_normalize(b_gene$log2FoldChange)) %>%
  dplyr::select(gene_id, allScore, log2FoldChange)


f_queen_score <- f_gene %>%
  mutate(allScore = min_max_normalize(f_gene$log2FoldChange)) %>%
  dplyr::select(gene_id, allScore, log2FoldChange)

f_worker_score <- f_gene %>%
  mutate(allScore = max_min_normalize(f_gene$log2FoldChange)) %>%
  dplyr::select(gene_id, allScore, log2FoldChange)


o_queen_score <- o_gene %>%
  mutate(allScore = min_max_normalize(o_gene$log2FoldChange)) %>%
  dplyr::select(gene_id, allScore, log2FoldChange)

o_worker_score <- o_gene %>%
  mutate(allScore = max_min_normalize(o_gene$log2FoldChange)) %>%
  dplyr::select(gene_id, allScore, log2FoldChange)


s_queen_score <- s_gene %>%
  mutate(allScore = min_max_normalize(s_gene$log2FoldChange)) %>%
  dplyr::select(gene_id, allScore, log2FoldChange)

s_worker_score <- s_gene %>%
  mutate(allScore = max_min_normalize(s_gene$log2FoldChange)) %>%
  dplyr::select(gene_id, allScore, log2FoldChange)

## Save results
write.csv(b_queen_score, file = "GO_data/brain_gene_queen_score.csv", row.names = FALSE)
write.csv(b_worker_score, file = "GO_data/brain_gene_worker_score.csv", row.names = FALSE)

write.csv(f_queen_score, file = "GO_data/fatbody_gene_queen_score.csv", row.names = FALSE)
write.csv(f_worker_score, file = "GO_data/fatbody_gene_worker_score.csv", row.names = FALSE)

write.csv(o_queen_score, file = "GO_data/ovary_gene_queen_score.csv", row.names = FALSE)
write.csv(o_worker_score, file = "GO_data/ovary_gene_worker_score.csv", row.names = FALSE)

write.csv(s_queen_score, file = "GO_data/spermatheca_gene_queen_score.csv", row.names = FALSE)
write.csv(s_worker_score, file = "GO_data/spermatheca_gene_worker_score.csv", row.names = FALSE)


