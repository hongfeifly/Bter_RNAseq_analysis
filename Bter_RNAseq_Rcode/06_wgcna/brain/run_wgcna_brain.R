#####
# author: "Joe Colgan", "Hongfei Xu"
# title: "**Weight gene co-expression network analysis of brain of
#   _Bombus terrestris_ workers and queens**"
# output:
#   pdf_document: default
#   html_document: default
#   fig_width: 4
#   fig_height: 4
#   fontsize: 20pt
#####

# Set the work directory
setwd("E:/rna_seq_analysis/results/Advanced_analysis/WGCNA/brain")
getwd()

## load the libraries
library(tidyverse)
library(WGCNA)
library(DESeq2)
library(igraph)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

## Load experimental data:
data <- readRDS(file = "../data/vsd_brain_caste.rds")
data_counts <- assay(data)

## Remove unexpected rows
data_counts <- data_counts %>%
  subset(!rownames(data_counts) %in% c("N_unmapped",
                                  "N_multimapping",
                                  "N_noFeature",
                                  "N_ambiguous"))

## We need a transposed dataframe consisting of gene per column and sample per row:
transposed_data <- t(data_counts)

## Run for experimental dataset:
transposed_gsg <- goodSamplesGenes(transposed_data, verbose = 3)
transposed_gsg$allOK

## Run for experimental samples:
exp_sampleTree = hclust(dist(transposed_data), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(exp_sampleTree, main = "Sample clustering to detect outliers",
     sub = "",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2)

## For the experimental dataset, there is an issue with sample BtQc_B4:
# Plot a line to show the cut
#abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(exp_sampleTree,
                     cutHeight = 60,
                     minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
transposed_data_filt = transposed_data[keepSamples, ]
nGenes = ncol(transposed_data_filt)
nSamples = nrow(transposed_data_filt)

## Read in sample information:
samples_information <- read.csv(file = "../data/brain_sample_information_WGCNA.csv",
                                row.names = 1)

# Save WGCNA basic data for use in subsequent parts
save(transposed_data_filt,
     samples_information,
     file = "results/brain_basic_data.RData")

## Run for experimental dataset:
# Re-cluster samples
exp_sampleTree2 = hclust(dist(transposed_data_filt),
                         method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
exp_traitColors = numbers2colors(samples_information, signed = FALSE);

# Plot the sample dendrogram and the colors underneath:
pdf(file = "results/sample_dendrogram_trait_heatmap_brain.pdf",
    height = 6,
    width = 15)
plotDendroAndColors(exp_sampleTree2,
                    exp_traitColors,
                    groupLabels = names(samples_information),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

saveRDS(object = exp_sampleTree2, file = "results/exp_sampleTree2_brain.rds")
saveRDS(object = exp_traitColors, file = "results/exp_traitColors_brain.rds")

## 2.b Step-by-step network construction and module detectio
## 2.b.1 Choosing the soft-thresholding power: analysis of network topology
## For the experimental set:
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12,
                        to = 20,
                        by = 2))
# Call the network topology analysis function
exp_sft = pickSoftThreshold(transposed_data_filt,
                            powerVector = powers,
                            verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(exp_sft$fitIndices[,1],
     -sign(exp_sft$fitIndices[,3])*exp_sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n",
     main = paste("Scale independence"));
text(exp_sft$fitIndices[,1],
     -sign(exp_sft$fitIndices[,3])*exp_sft$fitIndices[,2],
     labels = powers,
     cex = cex1,
     col = "red");
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.85,
       col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(exp_sft$fitIndices[,1],
     exp_sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(exp_sft$fitIndices[,1],
     exp_sft$fitIndices[,5],
     labels = powers,
     cex = cex1,
     col = "red")
abline(h = 100,
       col = "red")

exp_sft$powerEstimate
## exp_sft$fitIndices

## 2.b.2 Co-expression similarity and adjacency
## For experimental dataset:
softPower = 8;

################################################################################
Connectivity=softConnectivity(transposed_data_filt,corFnc = "cor", corOptions = "use ='p'",power=softPower,type = "signed")

## Plot connectivity, the slope should less than 0
pdf("results/scale-free connectivity.pdf");

scaleFreePlot(Connectivity,nBreaks = 10,truncated = FALSE,removeFirst = FALSE, main = "");

dev.off()
################################################################################

exp_adjacency = adjacency(transposed_data_filt,
                          power = softPower);

## 2.b.3 Topological Overlap Matrix (TOM)
## Experimental dataset:
# Turn adjacency into topological overlap
exp_TOM = TOMsimilarity(exp_adjacency);
diss_exp_TOM = 1-exp_TOM

## 2.b.4 Clustering using TOM
## Experimental data:
exp_geneTree = hclust(as.dist(diss_exp_TOM),
                      method = "average");

## Experimental set:
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
exp_dynamicMods = cutreeDynamic(dendro = exp_geneTree,
                                distM = diss_exp_TOM,
                                deepSplit = 2,
                                pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize);
table(exp_dynamicMods)

## Experimental data:
# Convert numeric lables into colors
exp_dynamicColors = labels2colors(exp_dynamicMods)
table(exp_dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(exp_geneTree,
                    exp_dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

## 2.b.5 Merging of modules whose expression profiles are very similar
## Experimental dataset:
# Calculate eigengenes
exp_MEList = moduleEigengenes(transposed_data_filt,
                              colors = exp_dynamicColors)
exp_MEs = exp_MEList$eigengenes
# Calculate dissimilarity of module eigengenes
exp_MEDiss = 1-cor(exp_MEs);
# Cluster module eigengenes
exp_METree = hclust(as.dist(exp_MEDiss),
                    method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(exp_METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")
abline(h = 0.25,
       col = "red")
## for the experimental dataset:
exp_MEDissThres = 0.25
# Plot the cut line into the dendrogram
# abline(h=exp_MEDissThres, col = "red")
# Call an automatic merging function
exp_merge = mergeCloseModules(transposed_data_filt,
                              exp_dynamicColors,
                              cutHeight = exp_MEDissThres,
                              verbose = 3)
# The merged module colors
exp_mergedColors = exp_merge$colors;
# Eigengenes of the new merged modules:
exp_mergedMEs = exp_merge$newMEs;

sizeGrWindow(12, 9)
pdf(file = "results/geneDendro-brain.pdf", wi = 9, he = 6)
plotDendroAndColors(exp_geneTree, cbind(exp_dynamicColors,
                                        exp_mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

dev.off()
## Experimental dataset:
# Rename to moduleColors
exp_moduleColors = exp_mergedColors
# Construct numerical labels corresponding to the colors
exp_colorOrder = c("grey", standardColors(50));
exp_moduleLabels = match(exp_moduleColors,
                         exp_colorOrder)-1;
exp_MEs = exp_mergedMEs;
dir.create(path = "results")
# Save module colors and labels for use in subsequent parts
save(exp_MEs,
     exp_moduleLabels,
     exp_moduleColors,
     exp_geneTree,
     file = "results/brain-02-networkConstruction-stepByStep.RData")

## 3 Relating modules to external clinical traits 3.a Quantifying module–trait associations
# Define numbers of genes and samples
exp_nGenes = ncol(transposed_data_filt);
exp_nSamples = nrow(transposed_data_filt);
# Recalculate MEs with color labels
exp_MEs0 = moduleEigengenes(transposed_data_filt,
                            exp_moduleColors)$eigengenes
exp_MEs = orderMEs(exp_MEs0)
exp_moduleTraitCor = cor(exp_MEs,
                         samples_information,
                         use = "p");
exp_moduleTraitPvalue = corPvalueStudent(exp_moduleTraitCor,
                                         exp_nSamples);
## Experimental dataset:
sizeGrWindow(10,6)
# Will display correlations and their p-values
exp_textMatrix = paste(signif(exp_moduleTraitCor, 2), "\n(",
                       signif(exp_moduleTraitPvalue, 1), ")", sep = "");
dim(exp_textMatrix) = dim(exp_moduleTraitCor)
# Display the correlation values within a heatmap plot:
pdf(file = "results/module-trait_relationship_heatmap_brain.pdf",
    height = 6,
    width = 8)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = exp_moduleTraitCor,
               xLabels = names(samples_information),
               yLabels = names(exp_MEs),
               ySymbols = names(exp_MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = exp_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# Save module colors and labels for use in subsequent parts
save(exp_MEs,
     exp_moduleTraitCor,
     samples_information,
     exp_textMatrix,
     file = "results/brain-03-input_for_labelled_heatmap.RData")

## 3.b Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable treatment containing the treatment column of samples information
caste = as.data.frame(samples_information$caste);
names(caste) = "caste"
# names (colors) of the modules
exp_modNames = substring(names(exp_MEs), 3)
exp_geneModuleMembership = as.data.frame(cor(transposed_data_filt,
                                             exp_MEs,
                                             use = "p"));
exp_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(exp_geneModuleMembership),
                                              exp_nSamples));

names(exp_geneModuleMembership) = paste("MM",
                                        exp_modNames,
                                        sep="");
names(exp_MMPvalue) = paste("p.MM",
                            exp_modNames,
                            sep = "");
exp_geneTraitSignificance = as.data.frame(cor(transposed_data_filt,
                                              caste,
                                              use = "p"));
exp_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(exp_geneTraitSignificance),
                                              exp_nSamples));
names(exp_geneTraitSignificance) = paste("GS.",
                                         names(caste),
                                         sep = "");
names(exp_GSPvalue) = paste("p.GS.",
                            names(caste),
                            sep = "")

## 3.c Intramodular analysis: identifying genes with high GS and MM
exp_module = "pink"
exp_column = match(exp_module,
                   exp_modNames);
exp_moduleGenes = exp_moduleColors == exp_module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(exp_geneModuleMembership[exp_moduleGenes, exp_column]),
                   abs(exp_geneTraitSignificance[exp_moduleGenes, 1]),
                   xlab = paste("Module Membership in", exp_module, "module"),
                   ylab = "Gene significance for caste",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2,
                   cex.lab = 1.2,
                   cex.axis = 1.2,
                   col = exp_module)

# Save module colors and labels for use in subsequent parts
save(exp_geneModuleMembership,
     exp_moduleGenes,
     exp_module,
     exp_column,
     exp_geneTraitSignificance,
     file = "results/brain-04-input_for_membership_scatterplot_pink.RData")

################################################################################
## Further analysis
## Load experimental data:
load(file = "results/brain_basic_data.RData")
load(file = "results/brain-02-networkConstruction-stepByStep.RData")
load(file = "results/brain-04-input_for_membership_scatterplot_pink.RData")

## add a column of gene id for MM and GS
exp_mm <- exp_geneModuleMembership %>%
  mutate(gene_id = rownames(exp_geneModuleMembership))
caste_sig <- exp_geneTraitSignificance %>%
  mutate(gene_id = rownames(exp_geneTraitSignificance))
colnames(caste_sig)[1] <- "GS_caste"
## combine GS and MM values together for each interested module
# and save results
exp_moduleGenes = exp_moduleColors == "pink"
b_pink_data <- exp_mm[exp_moduleGenes, c("gene_id", "MMpink")] %>%
  left_join(caste_sig, by = "gene_id")

exp_moduleGenes = exp_moduleColors == "white"
b_white_data <- exp_mm[exp_moduleGenes, c("gene_id", "MMwhite")] %>%
  left_join(caste_sig, by = "gene_id")

exp_moduleGenes = exp_moduleColors == "darkorange"
b_do_data <- exp_mm[exp_moduleGenes, c("gene_id", "MMdarkorange")] %>%
  left_join(caste_sig, by = "gene_id")

exp_moduleGenes = exp_moduleColors == "blue"
b_blue_data <- exp_mm[exp_moduleGenes, c("gene_id", "MMblue")] %>%
  left_join(caste_sig, by = "gene_id")

save(b_pink_data,
     b_white_data,
     b_do_data,
     b_blue_data,
     file = "results/brain_MM_GS_data.RData")
