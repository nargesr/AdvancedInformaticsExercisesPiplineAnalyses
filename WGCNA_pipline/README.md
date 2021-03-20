# Project

## Data description

We have 28 different samples which half of them are wild type and hal of the are [3xTgAD mouse](https://www.alzforum.org/research-models/3xtg). I used STAR and rsem to create RSEM file, then I extract TPM(Transcripts per million) as an input for our analysis.

After I removed the outlier based on [this figure](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/blob/main/WGCNA_pipline/sampleClusteringCleaning.pdf), I choose the power 14 to have scale free network [here](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/blob/main/WGCNA_pipline/summarypower.pdf) and find cluster base on expression and merge similar cluster [here](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/blob/main/WGCNA_pipline/geneDendro-3.pdf).

Now we have our cluster which means we can anaylse them! For analysing, firste I calculate the correlation and p-value between each culster and feature of data like genotype [here](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/blob/main/WGCNA_pipline/Module-traitRelationships.pdf). Then I did gene ontology by Enrich R and pathway analysis by reactome for some clusters which is intersting for us (like magenta and and paleturquoise which has significant correlation with genotype (3xTgADHO and 3xTgWT) and put it in 'GO_pathwayAnalyses' directory.

