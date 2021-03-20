# Project

## Data description

We have 28 different samples which half of them are wild type and hal of the are [3xTgAD mouse](https://www.alzforum.org/research-models/3xtg). all the data are extracted from hippocampuse in female mouse in three different time points (4 month, 12 month and 18 month).

## Prepare data

For aligning data I used STAR and for quantification I used RSEM, then I extract TPM(Transcripts per million) as an gene expression for our analysis.

## Analysis

After I checked my data to remove the outlier based on [this figure](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/blob/main/WGCNA_pipline/sampleClusteringCleaning.pdf), I choose the power 14 to have scale free network [here](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/blob/main/WGCNA_pipline/summarypower.pdf) and find cluster base on gene expression and merge similar cluster [here](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/blob/main/WGCNA_pipline/geneDendro-3.pdf). At the end, I have my finalized culters.

For analysing cluster, firste I calculate the correlation and p-value between each culster and feature of data like genotype [here](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/blob/main/WGCNA_pipline/Module-traitRelationships.pdf). Then I did gene ontology by Enrich R and pathway analysis by reactome for some clusters which is intersting for us (like magenta and and paleturquoise which has significant correlation with genotype (3xTgADHO and 3xTgWT) and put it in 'GO_pathwayAnalyses' directory.

For more infromation you can see the [WGCNA utorial](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)
