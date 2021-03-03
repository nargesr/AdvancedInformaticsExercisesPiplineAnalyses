#!/bin/bash
#SBATCH --job-name=countReadRNA      ## Name of the job.
#SBATCH -A ecoevo283         ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-20         ## number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs


module load subread/2.0.1

gtf="/data/class/ecoevo283/nargesr/ref/dmel-all-r6.13.gtf"
inpath="/data/class/ecoevo283/nargesr/RNAseq/"
file=$inpath"prefixes.txt" # prefix file
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`

# -p paired end, -T ?? num threads, -t Specify feature type(s) in a GTF annotation, -g Specify attribute type in GTF annotation, -Q minimum mapping quality score, -F format of the provided annotation file, -a gene annot, -o output file name

featureCounts -p -T 5 -t exon -g gene_id -Q 30 -F GTF -a $gtf -o ${inpath}counts/${prefix}.counts.txt ${inpath}mapped/${prefix}.sorted.bam
