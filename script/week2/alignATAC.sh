#!/bin/bash
#SBATCH --job-name=alignATAC ## Name of the job.
#SBATCH -A ecoevo283         ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-20         ## number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs
#SBATCH --output=alignATAC-%J.out ## output log file
#SBATCH --error=alignATAC-%J.err ## error log file


module load bwa/0.7.8
module load samtools/1.10
module load bcftools/1.10.2
module load python/3.8.0
module load gatk/4.1.9.0
module load picard-tools/1.87
module load java/1.8.0
module load hisat2/2.2.1

ls /data/class/ecoevo283/nargesr/ATACseq/data/*_R1.fq.gz | sed 's/\/data\/class\/ecoevo283\/nargesr\/ATACseq\/data\///; s/_R1.fq.gz//' > /data/class/ecoevo283/nargesr/ATACseq/prefixes.txt

inpath="/data/class/ecoevo283/nargesr/ATACseq/"
file=$inpath"prefixes.txt" # prefix file
ref="/data/class/ecoevo283/nargesr/ref/dmel-all-chromosome-r6.13.fasta"
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`

# alignments
bwa mem -t 2 -M $ref ${inpath}data/${prefix}_R1.fq.gz ${inpath}data/${prefix}_R2.fq.gz | samtools view -bS - > ${inpath}mapped/${prefix}.bam
samtools sort ${inpath}mapped/${prefix}.bam -o ${inpath}mapped/${prefix}.sort.bam
