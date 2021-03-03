#!/bin/bash
#SBATCH --job-name=callSNPDNA      ## Name of the job.
#SBATCH -A ecoevo283         ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-20         ## number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs


module load java/1.8.0
module load gatk/4.1.9.0
module load python/3.8.0

## create 10Mb regions
ref="/data/class/ecoevo283/nargesr/ref/dmel-all-chromosome-r6.13.fasta"
python3 fasta_generate_regions.py $ref 10000000 >/data/class/ecoevo283/nargesr/ref/my_regions_4Mb.txt

cd /data/class/ecoevo283/erebboah/DNAseq/gatk
interval=`head -n $SLURM_ARRAY_TASK_ID my_regions_4Mb.txt | tail -n 1`
/opt/apps/gatk/4.1.9.0/gatk GenotypeGVCFs -R $ref -V allsample.g.vcf.gz --intervals $interval -stand-call-conf 5 -O SNPbyregion/$interval.vcf.gz
## end second approach to calling SNPs