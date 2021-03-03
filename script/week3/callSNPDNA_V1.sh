#!/bin/bash
#SBATCH --job-name=callSNPDNA      ## Name of the job.
#SBATCH -A ecoevo283         ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-20         ## number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs


module load java/1.8.0
module load gatk/4.1.9.0

ref="/data/class/ecoevo283/nargesr/ref/dmel-all-chromosome-r6.13.fasta"
cd /data/class/ecoevo283/erebboah/DNAseq/gatk
/opt/apps/gatk/4.1.9.0/gatk GenotypeGVCFs -R $ref -V allsample.g.vcf.gz -stand-call-conf 5 -O result.vcf.gz
## end call SNPs


