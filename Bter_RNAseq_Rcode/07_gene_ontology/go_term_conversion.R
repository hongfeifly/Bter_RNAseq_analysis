#####
# author: "Joe Colgan" (joscolgan), "Hongfei Xu" (hongfeifly)
# title: "**go_term conversion**"
# output:
#   pdf_document: default
#   html_document: default
#   fig_width: 4
#   fig_height: 4
#   fontsize: 20pt
# descriptions: Due to the shortage of functional information for genes in
#   B. terrestris reference genome assembly, we obtained GO terms for
#   homologous genes in Drosophila melanogaster via Ensembl Metazoa BioMarts.
#####

## set directory
setwd("E:/rna_seq_analysis/results/Go")

## Check working directory:
getwd()

## load libraries
library(tidyverse)

## Load data into dataframe:
go_term_input <- read.table(file = "data/mart_export.txt",
                            stringsAsFactors = FALSE,
                            sep = "\t",
                            header = TRUE)

## Aggregate the list, which concatenates the second column:
mymerge <- function(x) {
  all_in_one <- paste(unlist(x), sep=",", collapse=",")
  split_term <- unlist(strsplit(all_in_one, split=","))
  return(paste(unique(split_term), sep=",", collapse=","))
}

## Run function:
aggregated_terms <- aggregate(go_term_input[-1], 
                              by = list(go_term_input$Gene.stable.ID),
                              mymerge)

write.table(aggregated_terms,
            file = "results/go_terms.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
