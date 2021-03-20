# Project

## Data description

We have 28 different samples which half of them are wild type and hal of the are [3xTgAD mouse](https://www.alzforum.org/research-models/3xtg). I used STAR and rsem to create RSEM file, then I extract TPM(Transcripts per million) as an input for our analysis.

After I removed the outlier based on below figure, I choose the power to have scale free network and find cluster base on expression.

<image src="sampleClusteringCleaning.pdf"/>
