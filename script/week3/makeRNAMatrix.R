library(data.table)
setwd("/data/class/ecoevo283/nargesr/RNAseq")

samples = read.table("/data/class/ecoevo283/nargesr/RNAseq/prefixes.txt",header=FALSE)
samples = samples$V1

thisMat=as.data.frame(read.table(paste0("/data/class/ecoevo283/nargesr/RNAseq/counts/",samples[1],".counts.txt"),sep="",head=T, stringsAsFactors = F))
colnames(thisMat) = c("GeneID","Chr","Start","End","Strand","Length","SampleID")
counts=thisMat[,c('GeneID','SampleID')]

for(i in 2:length(samples)){
  thisMat=as.data.frame(read.table(paste0("/data/class/ecoevo283/nargesr/RNAseq/counts/",samples[i],".counts.txt"),sep="",head=T, stringsAsFactors = F))
  colnames(thisMat) = c("GeneID","Chr","Start","End","Strand","Length","SampleID")
  thisMat=thisMat[na.omit(match(counts$GeneID,thisMat$GeneID)),]
  counts=as.data.frame(cbind(counts,thisMat[,"SampleID"]))
}

rownames(counts) = counts$GeneID
counts$GeneID = NULL
colnames(counts) = samples

write.table(counts, file="/data/class/ecoevo283/nargesr/RNAseq/counts/counts_100samples_flybaseIDs.tsv", quote=F)
