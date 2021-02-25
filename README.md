# Advanced Informatics Exercises Winter 2021

##  Week 1

i) Organize "the data"

I use pythonc and symlinks (`os.symlink()`) to link the data. you can find scripts in [`script/week1/`](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/tree/main/script/week1).

ii)  do some qc on one raw data file using fastqc

I ran fastqc on ADL06_1_1 raw data from the DNA-seq experiment. You can find the result [`fastqc/`](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/tree/main/fastqc).


## Week 2

Align the various datasets to the reference genome

you can find the related scripts in [`script/week2/`](https://github.com/nargesr/AdvancedInformaticsExercisesPiplineAnalyses/tree/main/script/week2).

For aligning DNA-seq and ATAC-seq, I used BWA v.0.7.0 and especificlly for DNA-seq, I added RG (read group which indexed with samtools) tags using Picard 1.87 (and Java 1.8.0) to use in GATK. The fastq reads convert to sam files by alignning, which are converted to bam files and sorted by using samtools 1.10.

For alignning RNA-seq I used HISAT2 v.2.2.1 and follow the same procedure
