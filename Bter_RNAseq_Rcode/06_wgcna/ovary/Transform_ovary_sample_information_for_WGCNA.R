#####
# author: "Joe Colgan" (joscolgan), "Hongfei Xu" (hongfeifly)
# title: "**Transform ovary sample informatin for
#   Weight gene co-expression network analysis**"
# output:
#   pdf_document: default
#   html_document: default
#   fig_width: 4
#   fig_height: 4
#   fontsize: 20pt
#####

# Set the work directory
setwd("E:/rna_seq_analysis/results/Advanced_analysis/WGCNA/ovary")
getwd()

## load the libraries
library(tidyverse)

## Read in sample information:
samples_infor <- read.table(file = "../data/ovary_sample_information_star.txt",
                                  header = FALSE)
samples_infor$V1 <- gsub(pattern = "ReadsPerGene.out.tab",
                               replacement = "",
                               samples_infor$V1)
## Remove tissue:
samples_infor$V3 <- NULL
row.names(samples_infor) <- samples_infor$V1
samples_infor$V1 <- NULL
## Update column names:
colnames(samples_infor) <- c("caste",
                             "stage",
                             "treatment")
## Remove outlier:
samples_infor <- subset(samples_infor,
                        row.names(samples_infor) %in% row.names(transposed_data_filt))
## Transform the sample information to number
samples_information <- samples_infor
samples_information$caste <- gsub(pattern = "Worker",
                                  replacement = 0,
                                  samples_infor$caste) %>%
  gsub(pattern = "Queen",
       replacement = 1,
       samples_infor$caste)
## stage
samples_information$stage <- gsub(pattern = "Stage_IV",
                                  replacement = 4,
                                  samples_infor$stage) %>%
  gsub(pattern = "Stage_II",
       replacement = 2,
       samples_infor$stage) %>%
  gsub(pattern = "Stage_I",
       replacement = 1,
       samples_infor)
## treatment
samples_information$treatment <- gsub(pattern = "Control",
       replacement = 0,
       samples_infor$treatment) %>%
  gsub(pattern = "Insemination_group",
       replacement = 1,
       samples_infor$treatment) %>%
  gsub(pattern = "Insemination_with_dilutent",
       replacement = 2,
       samples_infor$treatment) %>%
  gsub(pattern = "Insemination_with_nothing",
       replacement = 3,
       samples_infor$treatment)
samples_information$caste <- as.numeric(samples_information$caste)
samples_information$stage <- as.numeric(samples_information$stage)
samples_information$treatment <- as.numeric(samples_information$treatment)
## save data of sample information for WGCNA
write.csv(samples_information,
          file = "../data/ovary_sample_information_WGCNA.csv")
