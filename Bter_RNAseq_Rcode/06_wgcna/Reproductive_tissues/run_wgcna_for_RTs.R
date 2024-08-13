#####
# author: "Joe Colgan" (joscolgan), "Hongfei Xu" (hongfeifly)
# title: "**Weight gene co-expression network analysis of spermathecae of
#   _Bombus terrestris_ workers and queens**"
# output:
#   pdf_document: default
#   html_document: default
#   fig_width: 4
#   fig_height: 4
#   fontsize: 20pt
#####

# Set the work directory
setwd("E:/rna_seq_analysis/results/Advanced_analysis/WGCNA/spermatheca")
getwd()

## load the libraries
library(tidyverse)
library(WGCNA)
library(DESeq2)
library(igraph)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

## Load experimental data:
data <- readRDS(file = "../data/vsd_spermatheca_caste.rds")
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

## For the experimental dataset, there is an issue with sample 37:
# Plot a line to show the cut
#abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(exp_sampleTree,
                     cutHeight = 150,
                     minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
transposed_data_filt = transposed_data[keepSamples, ]
nGenes = ncol(transposed_data_filt)
nSamples = nrow(transposed_data_filt)

## Read in sample information:
samples_information <- read.csv(file = "../data/spermatheca_sample_information_WGCNA.csv",
                                row.names = 1)

# Save WGCNA basic data for use in subsequent parts
save(transposed_data_filt,
     samples_information,
     file = "results/spermatheca_basic_data.RData")

## Run for experimental dataset:
# Re-cluster samples
exp_sampleTree2 = hclust(dist(transposed_data_filt),
                         method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
exp_traitColors = numbers2colors(samples_information, signed = FALSE);

# Plot the sample dendrogram and the colors underneath:
pdf(file = "results/sample_dendrogram_trait_heatmap_RTs.pdf",
    height = 6,
    width = 8)
plotDendroAndColors(exp_sampleTree2,
                    exp_traitColors,
                    groupLabels = names(samples_information),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

saveRDS(object = exp_sampleTree2, file = "results/exp_sampleTree2_spermatheca.rds")
saveRDS(object = exp_traitColors, file = "results/exp_traitColors_spermatheca.rds")

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

## 2.b.2 Co-expression similarity and adjacency
## For experimental dataset:
softPower = 7;

################################################################################
Connectivity=softConnectivity(transposed_data_filt,corFnc = "cor", corOptions = "use ='p'",power=softPower,type = "signed")

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
pdf(file = "results/geneDendro-spermatheca.pdf", wi = 9, he = 6)
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
     file = "results/spermatheca-02-networkConstruction-stepByStep.RData")

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
pdf(file = "results/module-trait_relationship_heatmap_RTs.pdf",
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
     file = "results/spermatheca-03-input_for_labelled_heatmap.RData")

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
exp_module = "turquoise"
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
save(caste,
     exp_modNames,
     exp_geneModuleMembership,
     exp_MMPvalue,
     exp_geneTraitSignificance,
     exp_GSPvalue,
     exp_moduleGenes,
     exp_module,
     exp_column,
     file = "results/spermatheca-04-input_for_membership_scatterplot_turquoise.RData")

################################################################################
## Further analysis

## Load experimental data:
load(file = "results/spermatheca_basic_data.RData")
load(file = "results/spermatheca-02-networkConstruction-stepByStep.RData")
load(file = "results/spermatheca-04-input_for_membership_scatterplot_turquoise.RData")


######### Save data for plot correlation between MM and GS
## add a column of gene id for MM and GS
exp_mm <- exp_geneModuleMembership %>%
  mutate(gene_id = rownames(exp_geneModuleMembership))
caste_sig <- exp_geneTraitSignificance %>%
  mutate(gene_id = rownames(exp_geneTraitSignificance))
colnames(caste_sig)[1] <- "GS_caste"
## combine GS and MM values together for each interested module
# and save results
exp_moduleGenes = exp_moduleColors == "turquoise"
tur_data <- exp_mm[exp_moduleGenes, c("gene_id", "MMturquoise")] %>%
  left_join(caste_sig, by = "gene_id")

# exp_moduleGenes = exp_moduleColors == "brown"
# brown_data <- exp_mm[exp_moduleGenes, c("gene_id", "MMbrown")] %>%
#   left_join(caste_sig, by = "gene_id")
# 
# exp_moduleGenes = exp_moduleColors == "black"
# black_data <- exp_mm[exp_moduleGenes, c("gene_id", "MMblack")] %>%
#   left_join(caste_sig, by = "gene_id")
# 
# exp_moduleGenes = exp_moduleColors == "blue"
# blue_data <- exp_mm[exp_moduleGenes, c("gene_id", "MMblue")] %>%
#   left_join(caste_sig, by = "gene_id")
# 
# save(tur_data,
#      brown_data,
#      black_data,
#      blue_data,
#      file = "results/spermatheca_MM_GS_data.RData")


########### Identifying most important genes
## for one determined characteristic inside of the cluster
geneInfo0 = data.frame(EST = colnames(transposed_data_filt),
                       moduleColor = exp_moduleColors,
                       exp_geneTraitSignificance,
                       exp_GSPvalue)
## rank genes based on MM and GS
modOrder = order(-abs(cor(exp_MEs, caste, use = "p")))
for (mod in 1:ncol(exp_geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, exp_geneModuleMembership[, modOrder[mod]], 
                         exp_MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", exp_modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", exp_modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.caste))
geneInfo = geneInfo0[geneOrder, ]
## if you want to write the information in a csv file, just uncomment line below
#write.csv(geneInfo, file = "results/geneInfo.csv")

## extract important genes in different module based on MM and GS value
hub <- abs(exp_geneModuleMembership$MMturquoise) > 0.9 &
  abs(exp_geneTraitSignificance) > 0.9
table(hub)

turquoise_hub <- subset(tur_data, abs(tur_data$MMturquoise) > 0.8 &
                          abs(tur_data$GS_caste) > 0.2)

## extract the first one important gene in each module
HubGenes <- chooseTopHubInEachModule(transposed_data_filt, exp_moduleColors)


########### If bsps grouped in same module
## Import data of BSPs infor
bsp <- read.csv("../plot_MM_GS/data/chr_BSPs_information_with_cluster.csv")
bsp_list <- unique(sort(c(bsp$gene1,
                          bsp$gene2)))

## extract bsps' information
bspinfo <- subset(geneInfo, EST %in% bsp_list)

## identify modules that bsps grouped
table(bspinfo$moduleColor)

dif_list <- unique(sort(c(subset(bsp, biased == "dif_biased")$gene1,
                          subset(bsp, biased == "dif_biased")$gene2)))
dif_info <- subset(geneInfo, EST %in% dif_list)
table(dif_info$moduleColor)

queen_list <- unique(sort(c(subset(bsp, biased == "queen_biased")$gene1,
                            subset(bsp, biased == "queen_biased")$gene2)))
queen_info <- subset(geneInfo, EST %in% queen_list)
table(queen_info$moduleColor)

worker_list <- unique(sort(c(subset(bsp, biased == "worker_biased")$gene1,
                             subset(bsp, biased == "worker_biased")$gene2)))
worker_info <- subset(geneInfo, EST %in% worker_list)
table(worker_info$moduleColor)

qn_list <- unique(sort(c(subset(bsp, biased == "queen_non_biased")$gene1,
                            subset(bsp, biased == "queen_non_biased")$gene2)))
qn_info <- subset(geneInfo, EST %in% qn_list)
table(qn_info$moduleColor)

wn_list <- unique(sort(c(subset(bsp, biased == "worker_non_biased")$gene1,
                             subset(bsp, biased == "worker_non_biased")$gene2)))
wn_info <- subset(geneInfo, EST %in% wn_list)
table(wn_info$moduleColor)

## 
caste_bsp <- read.csv("chr_BSP_pairs_caste_biases_infor.csv")

qb_list <- unique(sort(c(subset(caste_bsp, biased_1 == "Q")$gene_1,
                         subset(caste_bsp, biased_2 == "Q")$gene_2)))
qb_info <- subset(geneInfo, EST %in% qb_list)
table(qb_info$moduleColor)

wb_list <- unique(sort(c(subset(caste_bsp, biased_1 == "W")$gene_1,
                         subset(caste_bsp, biased_2 == "W")$gene_2)))
wb_info <- subset(geneInfo, EST %in% wb_list)
table(wb_info$moduleColor)

bsp_add <- bsp %>%
  left_join(geneInfo[, c(1,2)], by = c("gene1" = "EST")) %>%
  left_join(geneInfo[, c(1,2)], by = c("gene2" = "EST"))
bsp_add$module <- ifelse(bsp_add$moduleColor.x == bsp_add$moduleColor.y,
                         bsp_add$moduleColor.x,
                         "dif")
## count the numbers of BSPs sets located on same/dif chr
bsp_count <- bsp_add %>%
  group_by(infor) %>%
  summarise(tur = sum(module == "turquoise"),
            dif = sum(module == "dif"),
            black = sum(module == "black"),
            green = sum(module == "green"),
            grey = sum(module == "grey"))

bsp_count2 <- tmp4 %>%
  group_by(infor) %>%
  summarise(tur = sum(module == "turquoise"),
            dif = sum(module == "dif"),
            black = sum(module == "black"),
            green = sum(module == "green"),
            grey = sum(module == "grey"))
## count BSPs sets number
length(unique(bsp_count2$infor))
# 173 sets

## count the numbers of BSPs sets located on same/dif chr
table(bsp_chr_count$dif_count)
# same_counts: 151; dif_counts: 22
