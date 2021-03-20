library(WGCNA)
library(RColorBrewer)
library(Hmisc)
library("ComplexHeatmap")
options(stringsAsFactors = FALSE);

###---------------------- 
# Read data
print("Read Data")
expressionList = read.csv('../../Data/expressionList_3xTgAD_hipp_female', header = TRUE);

## Prepare and clean data
#Remove rows with less than 1 TPM
expressionList = expressionList[expressionList[,ncol(expressionList)]>1,]

datExpr0 = as.data.frame(t(expressionList[,-c(1)]));
names(datExpr0) = expressionList$gene_id;
rownames(datExpr0) = names(expressionList)[-c(1)];

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenes(datExpr0, verbose = 3);
#if not okay 
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## Clustering
sampleTree = hclust(dist(datExpr0), method = "average");
# The user should change the dimensions if the window is too large or too small.
pdf(file = "sampleClusteringCleaning.pdf", width = 25, height = 6);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 50000, col = "red");
dev.off();

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)

datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#convert trancript ID to gene ID
datExpr = as.data.frame(t(datExpr))
for (i in c(1:dim(datExpr)[1])) {
  row.names(datExpr)[i] = strsplit(row.names(datExpr)[i], "\\.")[[1]][1]
}
datExpr = as.data.frame(t(datExpr))

collectGarbage();

save(datExpr, file = "Data/data_input.RData")

###---------------------- 
## Modules construction
print("WGCNA")
enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
pdf(file = "summarypower.pdf", width = 10, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.891,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off();

## Set Power
softPower = 14;
adjacency = adjacency(datExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
save(TOM, file = "Data/TOM.RData")

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf(file = "dendrogram.pdf", width = 18, height = 6);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off();


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
# Plot the dendrogram and colors underneath
pdf(file = "DynamicTreeCut.pdf", width = 18, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off();


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf(file = "eigenesgenes.pdf", width = 24, height = 5);
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.3
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off();
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
pdf(file = "geneDendro-3.pdf", wi = 18, he = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off();

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, dynamicColors, file = "Data/Data-networkConstruction.RData")

###---------------------- 
## ANALYSING
print("Analysing")
load(file = "Data/data_input.RData");
datTraits = read.csv('Data/datTraits_3xTgAD_hipp_female', sep = ',', stringsAsFactors = FALSE);
load(file = "Data/TOM.RData");
load(file = "Data/Data-networkConstruction.RData");
annot = read.csv('../../../Data/geneList', sep = '\t', stringsAsFactors = FALSE);

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels 
datME = moduleEigengenes(datExpr, moduleColors)$eigengenes
datME$MEgrey = NULL
MEs = orderMEs(datME)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

### For paper
moduleTraitCor = moduleTraitCor[c(12,3,4,10,7,2,6,13,5,15,11,9,14,8,1),]
moduleTraitPvalue = moduleTraitPvalue[c(12,3,4,10,7,2,6,13,5,15,11,9,14,8,1),]


## names
xlabels = c("4 month", "12 month", "18 month", "3xTgADHO", "3xTgADWT")
ylabels = c()
for (label in rownames(moduleTraitCor)) {
  ylabels = c(ylabels, substr(label, 3, 10000))
}
ylabels = capitalize(ylabels)

## Quantifying module-trait associations
pdf(file = "Module-traitRelationships.pdf", width = 10, height = 8);
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 8, 0.5, 2));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = xlabels,
               yLabels = ylabels,
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               cex.lab = 1,
               font.lab.x = 2,
               font.lab.y = 10,
               main = "")
dev.off();



## Plot the dendrogram

pdf(file = "EigengeneDendrogram.pdf", width = 6, height = 4);
par(cex = 1.0)
plotEigengeneNetworks(datME, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off();

pdf(file = "EigengeneDendrogramHeatmap.pdf", width = 6, height = 4);
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off();

## Network visuallization
# Select module probes
modules = names(table(moduleColors))
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
cyt2 = exportNetworkToCytoscape(modTOM,nodeFile = paste("Data/all-nodes50genes_revised", ".txt",sep="" ), weighted=TRUE, threshold = 0,
                                nodeNames = modProbes,nodeAttr = moduleColors[inModule])
all_nodes = read.delim("Data/all-nodes50genes_revised.txt",sep = "\t")
table(all_nodes$nodeAttr.nodesPresent...)


for (module in modules) {
  # Select module probes
  probes = names(datExpr)
  inModule = is.finite(match(moduleColors, module));
  modProbes = probes[inModule];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  
  cytB = exportNetworkToCytoscape(modTOM,nodeFile = paste("Data/", paste(module, "Node.txt",sep="" ),sep="" ), weighted=TRUE, threshold = 0,
                                  nodeNames = modProbes,nodeAttr = moduleColors[inModule])
  
}


## Diagnostics: displaying module heatmap and the eigengene
## Add color bar
month = rep(0, dim(datExpr)[1])
tissue = rep("Hippocampus", dim(datExpr)[1])
sex = rep("Male", dim(datExpr)[1])
mouseLine = rep("5XfAD;BL6", dim(datExpr)[1])
for (i in c(1:dim(datExpr)[1])) {
  tmp = strsplit(row.names(datExpr)[i], "_")[[1]]
  if (tmp[4] == "cortex"){
    tissue[i] = "Cortex"
  }
  if (tmp[2] == "Female"){
    sex[i] = "Female"
  }
  if (tmp[3] == "3xTgADWT"){
    mouseLine[i] = "3xTgADWT"
  }
  if (tmp[3] == "3xTgADHO"){
    mouseLine[i] = "3xTgADHO"
  }
  if (tmp[1] == "X4mon"){
    month[i] = "4 month"
  }
  if (tmp[1] == "X8mon"){
    month[i] = "8 month"
  }
  if (tmp[1] == "X12mon"){
    month[i] = "12 month"
  }
  if (tmp[1] == "X18mon"){
    month[i] = "18 month"
  }
}

## Plot module heatmap eigen gene
modules = names(table(moduleColors))
modules = modules[modules != "grey"]
for (which.module in modules) {
  # Recalculate MEs with color labels
  label = t(scale(datExpr[, moduleColors==which.module]))
  colnames(label) = paste0("sample_", seq(ncol(label)))
  ME=datME[, paste("ME",which.module, sep="")]
  
  for (i in c(1:length(row.names(label)))) {
    index = which(row.names(label)[i] == annot$Gene.ID)
    if (length(index) != 0) {
      if (annot$Gene.Symbol[index] %in% c("Cst7", "Gfap", "Apoe", "Trem2", "Thy1", "Gm4924", "Gas5", "Tyrobp", "Ctsz", "Ccl6", "Ctss")) {
        row.names(label)[i] = annot$Gene.Symbol[index]
      }
      else {
        row.names(label)[i] = ""
      }
    } else {
      row.names(label)[i] = ""
    }
  }
  
  pdf(file = paste0("ModuleHeatmapEigengene_", which.module, ".pdf"), width = 15, height = 8);
  column_ha = HeatmapAnnotation(Timepoint = month,
                                Tissue = tissue,
                                Sex = sex,
                                Genotype = mouseLine,
                                col = list(Timepoint = c("4 month" = "cyan", "12 month" = "deepskyblue", 
                                                         "8 month" = "blue", "18 month" = "darkblue"),
                                           Tissue = c("Hippocampus" = "deeppink", "Cortex" = "darkviolet"),
                                           Sex = c("Female" = "green", "Male" = "yellow"),
                                           Genotype = c("3xTgADHO" = "purple", "3xTgADWT" = "lightpink")),
                                eigengeneExp = anno_barplot(ME, baseline = 0, gp = gpar(fill = which.module)),
                                show_annotation_name = c(Timepoint = F, Tissue = F, Sex = F, Genotype = F, eigengeneExpression = T),
                                gap = unit(2, "points"), show_legend = TRUE, annotation_height = c(1,1,1,1,5), height = unit(6, "cm"))
  heatmap = Heatmap(label, column_title = which.module, cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE,
                    show_column_names = FALSE, top_annotation = column_ha, column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                    show_heatmap_legend = FALSE)
  print(heatmap)
  
  dev.off();
}



## convert gene ID to gene Name for each module
modules = names(table(moduleColors))
for (module in modules) {
  node = read.csv(paste0("Data/", module, "Node.txt"), sep = "\t")
  colnames(node) = c("Gene ID", "Gene symbol", "Module class")
  for (i in c(1:dim(node)[1])) {
    node[i,1] = strsplit(node[i,1], "\\.")[[1]][1]
    if (length(which(annot$Gene.ID == node[i,1])) != 0) {
      node[i,2] = annot$Gene.Symbol[which(annot$Gene.ID == node[i,1])]
    }
  }
  write.csv(node, paste0("Data/", module, "Node.txt"), row.names=F, quote = F)
}

## find size of each module
modules = names(table(moduleColors))
len = c()
for (module in modules) {
  node = read.csv(paste("Data/", paste(module, "Node.txt",sep="" ), sep="" ))
  len = c(len, dim(node)[1])
}

len = cbind(modules, len)
len = as.data.frame(len)
colnames(len) = c("moduleColor", "size")
write.csv(len, 'Data/moduleSize.txt', row.names = FALSE, quote = F)


pdf(file = "sizeModule.pdf");
barplot(as.integer(len$size),
        #names.arg = len$moduleColor,
        col = len$moduleColor,
        main = "Size of each module",
        ylab = "number of genes",
        xlab = "modules")
dev.off()


###---------
## paper
## create barplot
data = matrix(0, 6, 5)
colnames(data) = c("Time", "Sex", "Tissue", "Genotype", "counts")
data = as.data.frame(data)

data$Sex = "Female"
data$Tissue = "Hippocampus"

data$Time[1] = "4 month"
data$Genotype[1] = "3xTgADHO"
data$counts[1] = 5

data$Time[2] = "4 month"
data$Genotype[2] = "3xTgADWT"
data$counts[2] = 5

data$Time[3] = "12 month"
data$Genotype[3] = "3xTgADHO"
data$counts[3] = 5

data$Time[4] = "12 month"
data$Genotype[4] = "3xTgADWT"
data$counts[4] = 5

data$Time[5] = "18 month"
data$Genotype[5] = "3xTgADHO"
data$counts[5] = 9

data$Time[6] = "18 month"
data$Genotype[6] = "3xTgADWT"
data$counts[6] = 9


pdf("fig1A.pdf", width = 6, height = 2)

plot = ggplot(data, aes(x=Time, y=counts, fill=Genotype)) +
  geom_bar(stat="identity", position="dodge2") +
  scale_fill_manual(values = c("3xTgADHO" = "purple", "3xTgADWT" = "lightpink", 
                               "Female" = "green",
                               "4 month" = "cyan", "12 month" = "deepskyblue", "18 month" = "darkblue",
                               "Hippocampus" = "deeppink")) +
  coord_flip() +
  theme_classic() +
  geom_text(aes(label=counts), hjust=-0.5, color="black",
            position = position_dodge(0.9), size=5) +
  geom_point(aes(x=0, y=0, fill="Female")) +
  geom_point(aes(x=0, y=0, fill="Hippocampus")) +
  geom_point(aes(x=0, y=0, fill="4 month")) +
  geom_point(aes(x=0, y=0, fill="12 month")) +
  geom_point(aes(x=0, y=0, fill="18 month")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour="black", size = 16),
        axis.title=element_text(size = 18, face = "bold"),
        legend.text=element_text(colour="black", size = 16),
        legend.title=element_text(colour="black", size = 16))

print(plot)

dev.off()

### fig 1D
df = datExpr
index = c()
for (i in c(1:dim(df)[1])) {
  tmp = strsplit(row.names(df)[i], "_")[[1]]
  if (tmp[3] == "3xTgADWT"){
    index = c(index, 1)
  }
  if (tmp[3] == "3xTgADHO"){
    index = c(index , 2)
  }
}

df = rbind(df[index == 2, ], df[index == 1, ])

dfME = rbind(datME[index == 2, ], datME[index == 1, ])

## Add color bar
month = rep(0, dim(df)[1])
tissue = rep("Hippocampus", dim(df)[1])
sex = rep("Male", dim(df)[1])
mouseLine = rep("5XfAD;BL6", dim(df)[1])
for (i in c(1:dim(df)[1])) {
  tmp = strsplit(row.names(df)[i], "_")[[1]]
  if (tmp[4] == "cortex"){
    tissue[i] = "Cortex"
  }
  if (tmp[2] == "Female"){
    sex[i] = "Female"
  }
  if (tmp[3] == "3xTgADWT"){
    mouseLine[i] = "3xTgADWT"
  }
  if (tmp[3] == "3xTgADHO"){
    mouseLine[i] = "3xTgADHO"
  }
  if (tmp[1] == "X4mon"){
    month[i] = "4 month"
  }
  if (tmp[1] == "X8mon"){
    month[i] = "8 month"
  }
  if (tmp[1] == "X12mon"){
    month[i] = "12 month"
  }
  if (tmp[1] == "X18mon"){
    month[i] = "18 month"
  }
}

## Plot module heatmap eigen gene
modules = names(table(moduleColors))
modules = modules[modules != "grey"]
for (which.module in modules) {
  # Recalculate MEs with color labels
  label = t(scale(df[, moduleColors==which.module]))
  colnames(label) = paste0("sample_", seq(ncol(label)))
  ME=dfME[, paste("ME",which.module, sep="")]
  
  for (i in c(1:length(row.names(label)))) {
    index = which(row.names(label)[i] == annot$Gene.ID)
    if (length(index) != 0) {
      if (annot$Gene.Symbol[index] %in% c("Cst7", "Gfap", "Apoe", "Trem2", "Thy1", "Gm4924", "Gas5", "Tyrobp", "Ctsz", "Ccl6", "Ctss")) {
        row.names(label)[i] = annot$Gene.Symbol[index]
      }
      else {
        row.names(label)[i] = ""
      }
    } else {
      row.names(label)[i] = ""
    }
  }
  
  pdf(file = paste0("ModuleHeatmapEigengene_", which.module, "_reorder.pdf"), width = 15, height = 8);
  column_ha = HeatmapAnnotation(Tissue = tissue,
                                Sex = sex,
                                Genotype = mouseLine,
                                Timepoint = month,
                                col = list(Timepoint = c("4 month" = "cyan", "12 month" = "deepskyblue", 
                                                         "8 month" = "blue", "18 month" = "darkblue"),
                                           Tissue = c("Hippocampus" = "deeppink", "Cortex" = "darkviolet"),
                                           Sex = c("Female" = "green", "Male" = "yellow"),
                                           Genotype = c("3xTgADHO" = "purple", "3xTgADWT" = "lightpink")),
                                eigengeneExp = anno_barplot(ME, baseline = 0, gp = gpar(fill = which.module)),
                                show_annotation_name = c(Timepoint = F, Tissue = F, Sex = F, Genotype = F, eigengeneExpression = T),
                                gap = unit(2, "points"), show_legend = TRUE, annotation_height = c(1,1,1,1,5), height = unit(6, "cm"))
  heatmap = Heatmap(label, column_title = which.module, cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE,
                    show_column_names = FALSE, top_annotation = column_ha, column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                    show_heatmap_legend = FALSE)
  print(heatmap)
  
  dev.off();
}
