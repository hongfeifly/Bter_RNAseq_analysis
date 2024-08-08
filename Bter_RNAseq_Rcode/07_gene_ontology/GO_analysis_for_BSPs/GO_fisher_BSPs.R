#####
# author: "Joe Colgan", "Hongfei Xu"
# title: "**Gene Ontology term enrichment analysis of BSPs of _Bombus terrestris_**"
# output:
#   pdf_document: default
#   html_document: default
#   fig_width: 4
#   fig_height: 4
#   fontsize: 20pt
## Introduction
# The purpose of the present analysis is to examine similarities and differences
# in terms of gene expression associated with the BSPs.

# This script is for gene ontology (GO) enrichment analysis of
# differentially expressed genes to explore the biological processes
# and molecular functions associated with DE genes, using the R package 'TopGO'
# (https://www.bioconductor.org/packages/3.7/bioc/vignettes/topGO/inst/doc/topGO.pdf).
# Twp input files are required for the running of the test:  
#   1) A genelist file:  
#   The genelist file is a tab-delimited file containing two columns:  
#   Column 1: Locus (contains gene or transcript name of interest).  
#   Column 2: Rank value of interest (e.g. log fold changes).

# 2) A GO database file:  
#   The GO database file is a tab-delimited file containing two columns:  
#   Column 1: Locus (contains gene or transcript name of interest).  
#   Column 2: Comma separated GO terms (e.g. GO:0000001, GO:0000002, etc.).  

# This script will prepare data for GO analysis and
# create a 'TopGO object' from which enrichment tests can be performed
# to explore GO terms significantly enriched within the dataset. 
# This script outputs a results table of significantly enriched GO terms.
#####

## setup work directory
setwd("E:/rna_seq_analysis/results/Go/GO_Fisher_paralogues")

## load libraries
library(topGO)

## Step One: Input files, define objects for running topGO: 
## define go term
go_term <- "../GO_data/go_terms.txt"

## Read in GO annotations:
gene2go_mapping <- readMappings(file = go_term)

gene_universe <- names(gene2go_mapping)

## Import gene of interest data
bom_para <- read.csv("../GO_data/bombus_paralogues_cluster.csv")

gene_of_interest <- as.character(unique(sort(c(bom_para$gene1,
                                               bom_para$gene2))))

## extract gene for analysis
genelist <- factor(as.integer(gene_universe %in% gene_of_interest))
names(genelist) <- gene_universe

## Define node size:
node_size <- 20

## setup output directory
out_dir <- "bom_para_results_20"

if (file.exists(out_dir)) {
  stop("The output directory:", out_dir, ", already exists",
       "Let's avoid overwriting")
} else {
  dir.create(out_dir,
             recursive = TRUE)
}

################################################################################
## Steps Two and Three: Create topGO Object & run tests for GO term enrichment
for (go_category in c("BP")) {
  # STEP TWO
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    allGenes    = genelist,
                    gene2GO     = gene2go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = node_size) # Modify to reduce/increase stringency.
  # STEP THREE
  ## Calculate fisher exact test using 'weight01' algorithm:
  BP_weight_fisher <- runTest(object    = my_go_data,
                              algorithm = "weight01",
                              statistic = "fisher")
  ## Combine results from statistical tests:
  BP_weight_output <- GenTable(object       = my_go_data,
                               weight_fisher = BP_weight_fisher,
                               orderBy       = "weight_fisher",
                               topNodes      = length(score(BP_weight_fisher)))
  ## Correct ks test for multiple testing:
  BP_weight_output$weight_fisher <- as.numeric(BP_weight_output$weight_fisher)
  BP_weight_output$weight_fisher_adjusted <- p.adjust(p = BP_weight_output$weight_fisher,
                                                      method = c("BH"))
  ## Subset calls with significance higher than expected:
  BP_weight_output_sig <- subset(x      = BP_weight_output,
                                 subset = (Significant > Expected))
  
  BP_weight_output_sig$GO.ID <- gsub(pattern = ":",
                                     replacement = "_",
                                     BP_weight_output_sig$GO.ID)
  
  ## Write to output:
  write.table(x         = BP_weight_output_sig,
              file      = file.path(out_dir,
                                    paste(go_category,
                                          "sig.tsv",
                                          sep = "_")),
              row.names = FALSE,
              sep       = "\t",
              quote = FALSE)
}



################################################################################
## Steps Two and Three: Create topGO Object & run tests for GO term enrichment
for (go_category in c("MF")) {
  # STEP TWO
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    allGenes    = genelist,
                    gene2GO     = gene2go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = node_size) # Modify to reduce/increase stringency.
  # STEP THREE
  ## Calculate fisher exact test using 'weight01' algorithm:
  MF_weight_fisher <- runTest(object    = my_go_data,
                              algorithm = "weight01",
                              statistic = "fisher")
  ## Combine results from statistical tests:
  MF_weight_output <- GenTable(object       = my_go_data,
                               weight_fisher = MF_weight_fisher,
                               orderBy       = "weight_fisher",
                               topNodes      = length(score(MF_weight_fisher)))
  ## Correct ks test for multiple testing:
  MF_weight_output$weight_fisher <- as.numeric(MF_weight_output$weight_fisher)
  MF_weight_output$weight_fisher_adjusted <- p.adjust(p = MF_weight_output$weight_fisher,
                                                      method = c("BH"))
  ## Subset calls with significance higher than expected:
  MF_weight_output_sig <- subset(x      = MF_weight_output,
                                 subset = (Significant > Expected))
  
  MF_weight_output_sig$GO.ID <- gsub(pattern = ":",
                                     replacement = "_",
                                     MF_weight_output_sig$GO.ID)
  ## Write to output:
  write.table(x         = MF_weight_output_sig,
              file      = file.path(out_dir,
                                    paste(go_category,
                                          "sig.tsv",
                                          sep = "_")),
              row.names = FALSE,
              sep       = "\t",
              quote = FALSE)
}


################################################################################
## Steps Two and Three: Create topGO Object & run tests for GO term enrichment
for (go_category in c("CC")) {
  # STEP TWO
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    allGenes    = genelist,
                    gene2GO     = gene2go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = node_size) # Modify to reduce/increase stringency.
  # STEP THREE
  ## Calculate fisher exact test using 'weight01' algorithm:
  CC_weight_fisher <- runTest(object    = my_go_data,
                              algorithm = "weight01",
                              statistic = "fisher")
  ## Combine results from statistical tests:
  CC_weight_output <- GenTable(object       = my_go_data,
                               weight_fisher = CC_weight_fisher,
                               orderBy       = "weight_fisher",
                               topNodes      = length(score(CC_weight_fisher)))
  ## Correct ks test for multiple testing:
  CC_weight_output$weight_fisher <- as.numeric(CC_weight_output$weight_fisher)
  CC_weight_output$weight_fisher_adjusted <- p.adjust(p = CC_weight_output$weight_fisher,
                                                      method = c("BH"))
  ## Subset calls with significance higher than expected:
  CC_weight_output_sig <- subset(x      = CC_weight_output,
                                 subset = (Significant > Expected))
  
  CC_weight_output_sig$GO.ID <- gsub(pattern = ":",
                                     replacement = "_",
                                     CC_weight_output_sig$GO.ID)
  
  ## Write to output:
  write.table(x         = CC_weight_output_sig,
              file      = file.path(out_dir,
                                    paste(go_category,
                                          "sig.tsv",
                                          sep = "_")),
              row.names = FALSE,
              sep       = "\t",
              quote = FALSE)
}
