load(file = "../Data/Data-networkConstruction.RData");

moduleColors = names(table(moduleColors))
moduleColors = moduleColors[moduleColors != "grey"]

remove(geneTree, MEs, dynamicColors, moduleLabels)

library(ggplot2)
library(plotly)
library(viridis)
library(enrichR)
library(gridExtra)
library(grid)
library(ggpubr)
library(Hmisc)
options(stringsAsFactors = FALSE);
library(ComplexHeatmap)


## GO Enrichr
websiteLive = TRUE
dbs = listEnrichrDbs()
if (is.null(dbs)) websiteLive = FALSE


for (moduleColor in moduleColors) {
  ## Reactome
  pathway = read.csv(paste0("Reactome/", moduleColor, ".csv"))
  pathway = pathway[, c(2,3,6,7)]
  
  pathway = pathway[pathway$Entities.FDR < 0.1,]
  pathway = pathway[pathway$Entities.pValue < 0.05,]
  
  pathway$Entities.pValue = -1 * log2(pathway$Entities.pValue)
  pathway$X.Entities.found = as.integer(pathway$X.Entities.found)
  pathway = pathway[order(pathway$Entities.pValue, decreasing = TRUE),]
  
  
  ## Enrichr
  node = read.csv(paste0("../Data/", moduleColor, "Node.txt"), header = T, sep = ",")
  dbs = c("Mouse_Gene_Atlas", "GO_Biological_Process_2018")
  enriched = NULL
  if (websiteLive) {
    enriched = enrichr(node$Gene.symbol, dbs)
  }
  goTerm = enriched[["GO_Biological_Process_2018"]]
  goTerm = goTerm[, c(1,2,4)]
  
  goTerm = goTerm[goTerm$Adjusted.P.value < 0.05, ]
  if(dim(goTerm)[1] != 0){
    for (i in c(1:dim(goTerm)[1])) {
      goTerm$Overlap[i] = strsplit(goTerm$Overlap[i], "/")[[1]][1]
    }
  }
  goTerm$Overlap = as.numeric(goTerm$Overlap)
  
  goTerm$minusLog.Adjusted.P.value = -1 * log2(goTerm$Adjusted.P.value)
  goTerm = goTerm[order(goTerm$minusLog.Adjusted.P.value, decreasing = TRUE), ]
  
  
  if (dim(pathway)[1] != 0) {
    pathway_plot = ggplot(head(pathway, n=20), aes(x=Entities.pValue, y=Pathway.name, size=X.Entities.found, color=Entities.FDR)) +
      geom_point(alpha=0.5) +
      scale_size(range = c(5, 15), name="Number of genes") +
      scale_color_viridis(name = "FDR", discrete=FALSE) +
      ggtitle("Pathway analysis by Reactome") +
      xlab("-log2(p-value)") +
      ylab("") + 
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 1),
            axis.title=element_text(size = 12, face = "bold", vjust=10),
            axis.text=element_text(size = 10, face = "bold", color = "black"),
            legend.text = element_text(size = 10, color = "black"),
            legend.title = element_text(size = 10, color = "black"))
  }
  
  
  if(dim(goTerm)[1] != 0){
    go_plot = ggplot(head(goTerm, n=20), aes(x=minusLog.Adjusted.P.value, y=Term, size=Overlap, color=Adjusted.P.value)) +
      geom_point(alpha=0.5) +
      scale_size(range = c(5, 15), name="Number of genes") +
      scale_color_viridis(name = "Adjusted p-value", discrete=FALSE) +
      ggtitle("Gene ontology by EnrichR") +
      xlab("-log2(p-value)") +
      ylab("") + 
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 1),
            axis.title=element_text(size = 12, face = "bold", vjust=10),
            axis.text=element_text(size = 10, face = "bold", color = "black"),
            legend.text = element_text(size = 10, color = "black"),
            legend.title = element_text(size = 10, color = "black"))
  }
  
  if (dim(pathway)[1] != 0 & dim(goTerm)[1] != 0) {
    pdf(paste0(moduleColor, "Modules.pdf"), height = round(min(12, max(5, 3*dim(pathway)[1]/4, 3*dim(goTerm)[1]/4))), width = 25)
    
    figure = annotate_figure(grid.arrange(pathway_plot, go_plot, nrow=1), top = text_grob(paste0(moduleColor, " module with ", dim(node)[1], " genes"), color = "black", face = "bold", size = 16))
    
    dev.off();
    
  } else if (dim(pathway)[1] == 0  & dim(goTerm)[1] != 0) {
    pdf(paste0(moduleColor, "Modules.pdf"), height = round(min(12, max(5, 3*dim(goTerm)[1]/4))), width = 10)
    
    go_plot = go_plot +
      ggtitle(paste0(moduleColor, " module with ", dim(node)[1], " genes\nGene ontology by EnrichR"))
    
    print(go_plot)
    
    dev.off();
  } else if (dim(pathway)[1] != 0  & dim(goTerm)[1] == 0) {
    pdf(paste0(moduleColor, "Modules.pdf"), height = round(min(12, max(5, 3*dim(pathway)[1]/4))), width = 10)
    
    pathway_plot =  pathway_plot +
      ggtitle(paste0(moduleColor, " module with ", dim(node)[1], " genes\nPathway analysis by Reactome"))
    
    print(pathway_plot)
    
    dev.off();
  }
  
  remove(goTerm, figure, go_plot, pathway, pathway_plot)
}


#######
# paper
df = NULL
term = NULL
for (moduleColor in c("paleturquoise", "darkolivegreen", "magenta")) {
  node = read.csv(paste0("../Data/", moduleColor, "Node.txt"), header = T, sep = ",")
  dbs = c("Mouse_Gene_Atlas", "GO_Biological_Process_2018")
  enriched = NULL
  if (websiteLive) {
    enriched = enrichr(node$Gene.symbol, dbs)
  }
  goTerm = enriched[["GO_Biological_Process_2018"]]
  goTerm = goTerm[, c(1,2,4)]
  
  goTerm = goTerm[goTerm$Adjusted.P.value < 0.05, ]
  if(dim(goTerm)[1] != 0){
    for (i in c(1:dim(goTerm)[1])) {
      goTerm$Overlap[i] = strsplit(goTerm$Overlap[i], "/")[[1]][1]
    }
  }
  goTerm$Overlap = as.numeric(goTerm$Overlap)
  
  goTerm$minusLog.Adjusted.P.value = -1 * log2(goTerm$Adjusted.P.value)
  goTerm = goTerm[order(goTerm$minusLog.Adjusted.P.value, decreasing = TRUE), ]
  
  goTerm$Adjusted.P.value = NULL
  
  if (dim(goTerm)[1] != 0){
    for (i in c(1:dim(goTerm)[1])) {
      goTerm$Term[i] = strsplit(goTerm$Term[i], " \\(")[[1]][1]
    }
  }
  
  
  if (dim(goTerm)[1] != 0){
    goTerm$modules = capitalize(moduleColor)
    if (is.null(df)){
      df = head(goTerm, n=5)
      term = head(goTerm, n=5)$Term
    } else {
      df = rbind(head(goTerm, n=5), df)
      term = c(head(goTerm, n=5)$Term, term)
    }
  }
}

df$Term = factor(df$Term, levels = term)

pdf("GOPaper.pdf", height = 5, width = 10)
go_plot = ggplot(df, aes(x=minusLog.Adjusted.P.value, y=Term)) +
          geom_point(aes(size=Overlap, color=modules), alpha=100) +
          labs(color="Module names", size="Number of genes") +
          scale_colour_manual(values = unique(df$modules)) +  
          scale_size_continuous(breaks = c(4, 9, 14, 19), range = c(2, 10)) +
          xlab("-log2(p-value)") +
          ylab("") + 
          theme(plot.title = element_text(size = 14, face = "bold", hjust = 1),
                              axis.title=element_text(size = 16, vjust=10),
                               axis.text=element_text(size = 16, color = "black"),
                              legend.text = element_text(size = 16, color = "black"),
                              legend.title = element_text(size = 16, color = "black"))
print(go_plot)

dev.off();


### heatmap
for (moduleColor in moduleColors) {
  node = read.csv(paste0("../Data/", moduleColor, "Node.txt"), header = T, sep = ",")
  dbs = c("Mouse_Gene_Atlas", "GO_Biological_Process_2018")
  enriched = NULL
  if (websiteLive) {
    enriched = enrichr(node$Gene.symbol, dbs)
  }
  goTerm = enriched[["GO_Biological_Process_2018"]]
  goTerm = goTerm[, c(1,4,9)]
  goTerm = goTerm[goTerm$Adjusted.P.value < 0.05, ]
  
  
  if(dim(goTerm)[1] != 0){
    output = data.frame()
    count = 1
    for (i in c(1:dim(goTerm)[1])) {
      genes = strsplit(goTerm$Genes[i], ";")[[1]]
      output[count, 1] = strsplit(goTerm$Term[i], " \\(")[[1]][1]
      goTerm$Term[i] = strsplit(goTerm$Term[i], " \\(")[[1]][1]
      for (j in c(1:length(genes))) {
        output[count, genes[j]] = 1
      }
      count = count + 1
    }
    row.names(output) = output$V1
    output$V1 = NULL
    output = as.matrix(output)
    output[is.na(output)] = 0
    
    output = head(output, n = 20)
    #if (dim(output)[2] > 50) {
    #  output = output[, c(1:50)]
    #}
    if (dim(output)[1] > 1) {
      index = c()
      for (i in c(1:dim(output)[2])) {
        if (all(output[, i] == rep(0, dim(output)[1]))) {
          index = c(index, FALSE)
        } else {
          index = c(index, TRUE)
        }
      }
      output = output[, index]
    }
    
    pValue = c()
    for (i in row.names(output)) {
      pValue = c(pValue, goTerm$Adjusted.P.value[goTerm$Term == i])
    }
    
    library(circlize)
    color = colorRamp2(c(0, 1), c("white", "red"))
    
    annotation_row = data.frame(pValue = pValue)
    
    pdf(paste0("GO/", moduleColor, ".pdf"), width = max(10, round(dim(output)[2]/2)), height = max(2, round(dim(output)[1]/2)))
    plot = pheatmap(output, col = color, legend = FALSE, cluster_rows = FALSE, cluster_cols = FALSE,
                    annotation_row = annotation_row, main = paste0(moduleColor, " modules (", dim(node)[1], "genes)"))
    print(plot)
    dev.off()
  }
}



