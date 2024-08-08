###
### dN/dS estimation for Bombus terrestris during divergence node with Bombus impatiens
###
#####
# author: "Joe Colgan", "Hongfei Xu"
# title: "**_dN/dS estimation for Bombus terrestris_**"
#
## Introduction
# The purpose of the present analysis is to calculate dN/dS value for 
# Bombus terrestris during divergence node with Bombus impatiens.
# Estimations during other three nodes (Bombus poliars, Apis mellifera, Habropoda laboriosa)
# are run with same method.

# Set the work directory
setwd("rna_seq_analysis/dNdS_estimation")

## load the libraries
library(tidyverse)
library(orthologr)

################################################################################
# Part 1: calculate dNdS value for genes in B. terrestris

# get a dNdS table using:
# 1) reciprocal best hit for orthology inference (RBH)
# 2) Needleman-Wunsch for pairwise amino acid alignments
# 3) pal2nal for codon alignments
# 4) Comeron for dNdS estimation
# 5) single core processing 'comp_cores = 1'
dnds_bt <- dNdS(query_file      = "data/Bombus_terrestris.Bter_1.0.cds.fa",
                subject_file    = "data/Bombus_impatiens.cds.fa",
                ortho_detection = "RBH", 
                aa_aln_type     = "pairwise",
                aa_aln_tool     = "NW", 
                codon_aln_tool  = "pal2nal", 
                dnds_est.method = "Comeron", 
                comp_cores      = 1)


## re_colname for the dnds result
colnames(dnds_bt) <- colnames(colname)

## Import data: gene id list of B.terrestris for all transcripts
id_list <- read.table("data/bt_geneid_list.txt",
                      sep = "\t")
colnames(id_list) <- c("t_id", "gene_id")

## get gene id for dnds data
dnds_id <- dnds_bt %>%
  left_join(id_list, by = c("query_id" = "t_id"))

## remove rows which contain NA
dnds_select <- na.omit(dnds_id)            

## save results
write.csv(dnds_id,
          file = "dNdS_estimation_btvsbi.csv",
          row.names = F)
