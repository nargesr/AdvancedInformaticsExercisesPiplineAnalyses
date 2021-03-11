library(data.table)
options(stringsAsFactors = FALSE);

path = "/data/class/ecoevo283/nargesr/RNAseq/"

####-----
# Create fly_counts_2 table
samples = read.table(paste0(path, "prefixes.txt"), header = FALSE)
samples = samples$V1

counts = read.csv(paste0(path, "counts/", samples[1], ".counts.txt"), sep="", header = T)
colnames(coutns) = c("GeneID", "Chr", "Start", "End", "Strand", "Length", "SampleID")
counts = counts[, c('GeneID', 'SampleID')]

for(i in 2:length(samples)){
  countsSample = read.csv(paste0(path, "counts/", samples[i], ".counts.txt"), sep="", head=T)
  colnames(countsSample) = c("GeneID", "Chr", "Start", "End", "Strand", "Length", "SampleID")
  countsSample = countsSample[na.omit(match(counts$GeneID, countsSample$GeneID)),]
  counts = as.data.frame(cbind(counts, countsSample[, "SampleID"]))
}

rownames(counts) = counts$GeneID
counts$GeneID = NULL
colnames(counts) = samples

write.table(counts, file=paste0(path, "/counts/fly_counts_2.tsv"), quote=F)


####-----
# create DEseq2 object & run DEseq

library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(EnhancedVolcano)

countData = read.table(paste0(path, "/counts/fly_counts_2.tsv"))
sampleInfo = read.table(paste0(path, "RNAseq384_SampleCoding.txt"), header=T)
sampleInfo_100 = sampleInfo[paste0("x",sampleInfo$FullSampleName) %in% colnames(countdata),]
sampleInfo_100 = sampleInfo_100[match(colnames(countdata), paste0("x",sampleInfo_100$FullSampleName)),]
dds = DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo_100, design=~TissueCode)
# filter lowly expressed genes
keep = rowSums(counts(dds)) >= 10
dds = dds[keep, ]
dds = DESeq(dds)
res = results(dds)

rld = rlog(dds)
mydata = assay(rld)

save(rld, file=paste0(path, "analysis/rld.rda"))
save(dds, file=paste0(path, "analysis/dds.rda"))


####----
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(EnhancedVolcano)

load(paste0(path, "analysis/rld.rda"))
load(paste0(path, "analysis/dds.rda"))
res = results( dds )

path = "/data/class/ecoevo283/erebboah/RNAseq/figures/"

pdf(file=paste0(path, "plotMA.pdf"), width = 5, height = 5) 
plotMA(res, ylim = c(-1, 1), colNonSig = "gray60", colSig = "blue")
dev.off()

pdf(file=paste0(path, "dispEst.pdf"), width = 5, height = 5) 
plotDispEsts(dds)
dev.off()

pdf(file=paste0(path, "hist.pdf"), width = 5, height = 5) 
hist( res$pvalue, breaks=20, col="grey" )
dev.off()

sampleDists = dist(t(assay(rld)))
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = rld$TissueCode
colnames(sampleDistMatrix) = NULL
colours = colorRampPalette( rev(brewer.pal(9, "Blues")))(255)

pdf(file=paste0(path, "heatmap.pdf"), width = 5, height = 5) 
heatmap.2( sampleDistMatrix, trace="none", col=colours)
dev.off()

pdf(file=paste0(path, "pca.pdf"), width = 5, height = 5) 
print(plotPCA(rld, intgroup = "TissueCode"))
dev.off()

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
pdf(file=paste0(path, "heatmapTopvargenes.pdf"), width =5, height = 5) 
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

res <- results(dds, contrast = c('TissueCode','P','B'))
res <- lfcShrink(dds, contrast = c('TissueCode','P','B'), res=res, type = 'normal')

pdf(file=paste0(path, "enhancedVolcano.pdf"), width =10, height = 10)     
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')
dev.off()

