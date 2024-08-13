#####
# author: Hongfei Xu (hongfeifly)
# title: "**gene extraction for tissue specificity calculation**"
# output:
#   pdf_document: default
#   html_document: default
#   fig_width: 4
#   fig_height: 4
#   fontsize: 20pt
# descriptions: The data was gained from STAR output and normalized by Deseq2,
#   then we calculated a tpm value for each gene, this analysis was ran for
#   removing the lowly expressed genes which under a tpm threshold.
#####

## Set work directory
setwd("E:/rna_seq_analysis/results/Tissue_specificity")

## Check working directory:
getwd()

## Load libraries:
library(tidyverse)
library(parallel)
library(GenomicFeatures)

## import the gtf file
txdb <- makeTxDbFromGFF("data/Bombus_terrestris.Bter_1.0.55.gtf",
                        format = "gtf")

## extact exons of every gene
exons_gene <- exonsBy(txdb, by = "gene")

## set up cl number, we will use half of the cores in our computer
cl <- makeCluster(0.50 * detectCores())

## calculate the length of exons, remove the redundant parts by "reduce",
# collect the length by "width", and get the final length.
exons_lens <- parLapply(cl, exons_gene, function(x){sum(width(reduce(x)))})
exons_lens[1:10]

## transfer the list to the dataframe
geneid_efflen <- data.frame(gene_id = names(exons_lens),
                            efflen = as.numeric(exons_lens))

### Import the gene counts from STAR.
counts_all <- read.csv("data/normalized_counts_all.csv")
colnames(counts_all)[1] <- "gene_id"

## create a function for calculate TPM value by counts
countToTpm <- function(counts, effLen)
{ 
  rate <- log(counts) - log(effLen) 
  denom <- log(sum(exp(rate))) 
  exp(rate - denom + log(1e6)) 
} 

## merge the efflen to counts file
counts <- counts_all %>%
  subset(!gene_id %in% c("N_unmapped",
                         "N_multimapping",
                         "N_noFeature",
                         "N_ambiguous")) %>%
  left_join(geneid_efflen, by = "gene_id")

## calculate the TPM value for each gene
tpm <- with(counts, countToTpm(counts = counts[, 2:257],
                               effLen = counts$efflen))

## add two columns in tpm_df: gene_id and tpm_mean
tpm$tpm_mean <- rowMeans(tpm)
tpm$gene_id <- counts$gene_id

## check the value of tpm_mean
quantile(log2(tpm$tpm_mean),
         0.10)

ggplot(data = tpm,
       aes(x = log2(tpm_mean))) +
  geom_density()

## filter the data, remove the first 10% tpm value
tpm_keep <- tpm[, c("gene_id", "tpm_mean")] %>%
  subset(log2(tpm_mean) > -8.620517) %>%
  left_join(geneid_efflen, by = "gene_id")

write.csv(tpm_keep, file = "data/tpm_keep_star.csv")

lintr::lint("Calculate_TPM_for_STAR.R")
