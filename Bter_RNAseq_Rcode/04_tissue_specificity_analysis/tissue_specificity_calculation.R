#####
# author: "Hongfei Xu" (hongfeifly)
# title: "**tissue specificity calculation**"
# output:
#   pdf_document: default
#   html_document: default
#   fig_width: 4
#   fig_height: 4
#   fontsize: 20pt
# descriptions: We extracted genes based on tpm value,
#   and calculated the tau value (tissue specificity) using tispec.
#####

setwd("rna_seq_analysis/Tissue_specificity")

getwd()

## load the libraries
library(tidyverse)
library(tispec)
library(ggpubr)
library(eulerr)

## Import the gene counts mean of all and control samples.
cm_all <- read.csv("data/normalized_counts_mean_all.csv", row.names = 1)

## import tpm value to evaluate imported gene counts
tpm_keep <- read.csv("data/tpm_keep_star.csv", row.names = 1)

## get the counts of the genes which in tpm_keep(the lowest 10% were removed) 
cm_all_keep <- cm_all %>%
  subset(rownames(cm_all) %in% tpm_keep$gene_id)

# Separate several different groups
meanExp_all <- cm_all_keep[, c("brain_mean", "fatbody_mean",
                               "ovary_mean", "spermatheca_mean")]

# check input data
checkInput(cm_all_keep)

# Normalise the data
## A pseudo count is added to log2 transform the data, and a threshold of log2(0) is set.
# This means that all genes with input values below 1.0 will be defined as non-expressed.
# All genes not expressed in any tissue are then removed from the data.
log2Exp_all <- log2Tran(meanExp_all)
head(log2Exp_all)
qnExp_all <- quantNorm(log2Exp_all)
head(qnExp_all)

### Calculate Tissue Specificity
tauExp_all <- calcTau(qnExp_all)
head(tauExp_all)

## save the results: tau value and selected gene counts mean
write.csv(tauExp_all, file = "results/tau_all_star.csv")

lintr::lint("tissue_specificity_calculation.R")
