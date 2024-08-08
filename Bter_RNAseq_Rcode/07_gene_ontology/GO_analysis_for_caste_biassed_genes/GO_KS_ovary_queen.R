#####
# author: "Joe Colgan", "Hongfei Xu"
# title: "**Gene Ontology term enrichment analysis of queen-biassed genes
# in ovary of _Bombus terrestris_**"
# output:
#   pdf_document: default
#   html_document: default
#   fig_width: 4
#   fig_height: 4
#   fontsize: 20pt
## Introduction
# The purpose of the present analysis is to examine similarities and differences
# in terms of gene expression associated with the queen-biassed gene in ovary.

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
setwd("E:/rna_seq_analysis/results/Go/GO_KS_casteDEGs")

## load libraries
library(topGO)

## Step One: Input files, define objects for running topGO: 

## Read in information
oq_score <- read.csv(file = "../GO_data/ovary_gene_queen_score.csv")

## define go term
go_term <- "../GO_data/go_terms.txt"

## Read in GO annotations:
gene2go_mapping <- readMappings(file = go_term)

## Convert into topgo's genelist format:
topgo_genelist <- oq_score$allScore
names(topgo_genelist) <- oq_score$gene_id

## Get top 5% genes
quantile(topgo_genelist, 0.05)

topDiffGenes <- function(allScore){
  return(allScore <= 0.2133915)
}

sum(topDiffGenes(topgo_genelist))
# 579 genes here, oq_bias(187)

## Define node size:
node_size <- 20

## setup output directory
out_dir <- "ovary_queen_results_20"

if (file.exists(out_dir)) {
  stop("The output directory:", out_dir, ", already exists",
       "Let's avoid overwriting")
} else {
  dir.create(out_dir,
             recursive = TRUE)
}

################################################################################
## go_cotegory = "BP": Create topGO Object & run tests for GO term enrichment

## Steps Two and Three
for (go_category in c("BP")) {
  # STEP TWO
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    geneSel     = topDiffGenes,
                    allGenes    = topgo_genelist,
                    gene2GO     = gene2go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = node_size) # Modify to reduce/increase stringency.
  # STEP THREE
  ## Calculate ks test using 'weight01' algorithm:
  BP_weight_ks <- runTest(object = my_go_data,
                          algorithm = "weight01",
                          statistic = "ks")
  ## Combine results from statistical tests:
  BP_weight_output <- GenTable(object        = my_go_data,
                               weight_ks     = BP_weight_ks,
                               orderBy       = "weight_ks",
                               topNodes      = length(score(BP_weight_ks)))
  ## Subset calls with significance higher than expected:
  BP_weight_output_sig <- subset(x      = BP_weight_output,
                                 subset = (Significant > Expected) &
                                   (weight_ks < 0.05))
  ## Print to console one of the GO terms of interest to check the distribution of that GO term across ranked genes:
  print(showGroupDensity(object  = my_go_data,
                         whichGO = head(BP_weight_output_sig$GO.ID,
                                        n = 1),
                         ranks   = TRUE,
                         rm.one  = FALSE))
  
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
## go_cotegory = "MF": Create topGO Object & run tests for GO term enrichment

## Steps Two and Three
for (go_category in c("MF")) {
  # STEP TWO
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    geneSel     = topDiffGenes,
                    allGenes    = topgo_genelist,
                    gene2GO     = gene2go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = node_size) # Modify to reduce/increase stringency.
  # STEP THREE
  ## Calculate ks test using 'weight01' algorithm:
  MF_weight_ks <- runTest(object = my_go_data,
                          algorithm = "weight01",
                          statistic = "ks")
  ## Combine results from statistical tests:
  MF_weight_output <- GenTable(object       = my_go_data,
                               weight_ks     = MF_weight_ks,
                               orderBy       = "weight_ks",
                               topNodes      = length(score(MF_weight_ks)))
  ## Subset calls with significance higher than expected:
  MF_weight_output_sig <- subset(x      = MF_weight_output,
                                 subset = (Significant > Expected) &
                                   (weight_ks < 0.05))
  ## Print to console one of the GO terms of interest to check the distribution of that GO term across ranked genes:
  print(showGroupDensity(object  = my_go_data,
                         whichGO = head(MF_weight_output_sig$GO.ID,
                                        n = 1),
                         ranks   = TRUE,
                         rm.one  = FALSE))
  
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
## go_cotegory = "CC": Create topGO Object & run tests for GO term enrichment

## Steps Two and Three
for (go_category in c("CC")) {
  # STEP TWO
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    geneSel     = topDiffGenes,
                    allGenes    = topgo_genelist,
                    gene2GO     = gene2go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = node_size) # Modify to reduce/increase stringency.
  # STEP THREE
  ## Calculate ks test using 'weight01' algorithm:
  CC_weight_ks <- runTest(object = my_go_data,
                          algorithm = "weight01",
                          statistic = "ks")
  ## Combine results from statistical tests:
  CC_weight_output <- GenTable(object       = my_go_data,
                               weight_ks     = CC_weight_ks,
                               orderBy       = "weight_ks",
                               topNodes      = length(score(CC_weight_ks)))
  ## Subset calls with significance higher than expected:
  CC_weight_output_sig <- subset(x      = CC_weight_output,
                                 subset = (Significant > Expected) &
                                   (weight_ks < 0.05))
  ## Print to console one of the GO terms of interest to check the distribution of that GO term across ranked genes:
  print(showGroupDensity(object  = my_go_data,
                         whichGO = head(CC_weight_output_sig$GO.ID,
                                        n = 1),
                         ranks   = TRUE,
                         rm.one  = FALSE))
  
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
