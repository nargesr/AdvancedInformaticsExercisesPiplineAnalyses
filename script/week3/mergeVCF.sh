#!/bin/bash
#SBATCH --job-name=mergeVCF      ## Name of the job.
#SBATCH -A ecoevo283         ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-20         ## number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs


module load java/1.8.0
module load gatk/4.1.9.0

ref="/data/class/ecoevo283/nargesr/ref/dmel-all-chromosome-r6.13.fasta"

mkdir /data/class/ecoevo283/nargesr/DNAseq/gatk
cd /data/class/ecoevo283/erebboah/DNAseq/gatk

/opt/apps/gatk/4.1.9.0/gatk CombineGVCFs -R $ref $(printf -- '-V %s ' DNAseq2out/*.g.vcf.gz) -O allsample.g.vcf.gz
## and merge gVCFs