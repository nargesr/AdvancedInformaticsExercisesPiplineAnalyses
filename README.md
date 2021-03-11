# Advanced Informatics Exercises Winter 2021

##  Week 1

*Goal one*: Organize "the data"

I use python and symlinks (`os.symlink()`) to link the data. you can find scripts in [`script/week1/`](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/tree/main/script/week1).

The scripts create files with the new names for the samples are in the directories `DNAseq`, `ATACseq`, or `RNAseq` with the symlinks in sub-directories named `data`:
```
DNAseq/
    dna_samples.txt
    data/
        ADL06_1_1.fq.gz
        ADL06_1_2.fq.gz
        ADL06_2_1.fq.gz
        ADL06_2_2.fq.gz
        ADL06_3_1.fq.gz
        ADL06_3_2.fq.gz
        ...
ATACseq/
    atac_samples.txt
    data/
        P004_R1.fq.gz
        P004_R2.fq.gz 
        P005_R1.fq.gz
        P005_R2.fq.gz
        P006_R1.fq.gz
        P006_R2.fq.gz
        ...
RNAseq/
    rna_samples.txt
    data/
        x21001B0_R1.fq.gz
        x21001B0_R2.fq.gz
        x21001E0_R1.fq.gz
        x21001E0_R2.fq.gz
        x21001H0_R1.fq.gz
        x21001H0_R2.fq.gz
        ...
```

*Goal two*: Do some qc on one raw data file using fastqc

I ran fastqc on ADL06_1_1 raw data from the DNA-seq experiment. You can find the result [`fastqc/`](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/tree/main/fastqc).


## Week 2

*Goal*: Align the various datasets to the reference genome.

I made "prefix" files at the begining of each bash scripts with the provided code from Dr. Long's notes.

The prefix files are in the `DNAseq`, `ATACseq`, and `RNAseq` folders.

The bash scripts that made use of the prefix files are called `alignDNA.sh`, `alignATAC.sh`, and `alignRNA.sh`. For aligning DNA-seq and ATAC-seq, I used BWA v.0.7.0 and especificlly for DNA-seq, I added RG (read group which indexed with samtools) tags using Picard 1.87 (and Java 1.8.0) to use in GATK. The fastq reads convert to sam files by alignning, which are converted to bam files and sorted by using samtools 1.10.

For alignning RNA-seq I used HISAT2 v.2.2.1 and follow the same procedure

The output is in each subfolder named `mapped`:
```
DNAseq/
    dna_samples.txt
    prefixes.txt
    data/
    mapped/
        ADL06_1.RG.bam
        ADL06_1.RG.bam.bai
        ADL06_1.sort.bam
        ADL06_2.RG.bam
        ADL06_2.RG.bam.bai
        ADL06_2.sort.bam
        ADL06_3.RG.bam
        ADL06_3.RG.bam.bai
        ADL06_3.sort.bam
        ...
ATACseq/
    atac_samples.txt
    prefixes.txt
    data/
    mapped/
        P004.sort.bam
        P005.sort.bam
        P006.sort.bam
        ...
RNAseq/
    rna_samples.txt
    prefixes.txt
    prefixes_random.txt
    data/
    mapped/
        x21001B0.sorted.bam
        x21001B0.sorted.bam.bai
        x21001B0.sorted.bam
        x21001H0.sorted.bam.bai
        x21001H0.sorted.bam
        x21001H0.sorted.bam.bai
        ...
```

you can find the related scripts in [`script/week2/`](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/tree/main/script/week2).


##  Week 3

*Goal one*: continue grinding through the results of our last week pipeline which is work through the DNA-seq analysis pipeline to call SNPs

*Goal two*: Start analysing RNA-seq by generating counts per gene for RNA-seq data using the subread package.

## DNA-seq GATK pipeline
1. Create a file call `prefix2.txt` in the `DNAseq` folder with the sample names like this:
```
ADL06
ADL09
ADL10
ADL14
```

2. Merge within samples and call SNPs by running `mergeSNPDNA.sh` and put the output in `DNAseq/gatk`:
```
DNAseq/
    dna_samples.txt
    prefixes.txt
    prefixes2.txt
    data/
    mapped/
    gatk/
        ADL06.dedup.bam
        ADL06.dedup.bam.bai
        ADL06.dedup.bam.sbi
        ADL06.dedup.metrics.txt
        ADL06.g.vcf.gz
        ADL06.g.vcf.gz.tbi
        ...
```

3. Merge by sample vcfs into a single vcf `mergeVCF.sh` and put the output in `DNAseq/gatk/allsample.g.vcf.gz` and index in `DNAseq/gatk/allsample.g.vcf.gz.tbi`.

4. Call SNPs by running `callSNPDNA_V1.sh` and put the output in `DNAseq/gatk/result.vcf.gz` and index in `DNAseq/gatk/result.vcf.gz.tbi`

Also there is second approach to calling SNPs that you can run in parallel by running `callSNPDNA_V2.sh`. In this bash script, I first run `fasta_generate_regions.py` to create 10Mb regions and then use the `--intervals` option in GATK with the `my_regions_10Mb.txt` file to output a `.vcf.gz` and index for each 10Mb region (1,881 regions total).


The output is in `DNAseq/gatk/SNPbyregion`:
```
DNAseq/
    dna_samples.txt
    prefixes.txt
    prefixes2.txt
    data/
    mapped/
    gatk/
       SNPbyregion/
            ADL06.dedup.bam.bai
            ADL06.dedup.bam.sbi
            ADL06.dedup.metrics.txt
            ADL06.g.vcf.gz
            ADL06.g.vcf.gz.tbi
            ...
```


## RNA-seq 
I follow the same aporaoch I used for DNA-seq to count the number of reads per gene from the 100 sorted `BAM` files using `featureCounts` from the subreads package and it to in `countReadRNA.sh`.

The output is in `RNAseq/counts`:
```
RNAseq/
    rna_samples.txt
    prefixes.txt
    prefixes_random.txt
    data/
    mapped/
    counts/
        x21001B0.counts.txt
        x21001B0.counts.txt.summary
        x21001H0.counts.txt
        x21001H0.counts.txt.summary
        x21001P0.counts.txt
        x21001P0.counts.txt.summary
        ...
```

you can find the related scripts in [`script/week3/`](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/tree/main/script/week3).


##  Week 4

*Goal one*: Creat a "featureCounts" object called "fly_counts_2.txt"  -- rows are genes, columns are samplenames for RNA-seq (I forgot to do it last week so have to do it for next step)

*Goal two*: Analyze RNAseq

*Goal three*: Visualize ATACseq peaks
