#!/bin/bash
#SBATCH --job-name=mergeSNPDNA      ## Name of the job.
#SBATCH -A ecoevo283         ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-20         ## number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs


module load bwa/0.7.17
module load samtools/1.10
module load bcftools/1.10.2
module load python/3.8.0
module load java/1.8.0
module load gatk/4.1.9.0
module load picard-tools/1.87
module load bamtools/2.5.1        # bamtools merge is useful
module load freebayes/1.3.2       # fasta_generate_regions.py is useful
module load vcftools/0.1.16

inpath="/data/class/ecoevo283/nargesr/DNAseq/"
file=$inpath"prefixes2.txt" # prefix file
ref="/data/class/ecoevo283/nargesr/ref/dmel-all-chromosome-r6.13.fasta"
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`

# merge lanes within samples
cd /data/class/ecoevo283/nargesr/DNAseq/mapped
java -jar /opt/apps/picard-tools/1.87/MergeSamFiles.jar $(printf 'I=%s ' ${prefix}*.RG.bam) O=${prefix}.merge.bam
/opt/apps/gatk/4.1.9.0/gatk MarkDuplicatesSpark -I ${prefix}.merge.bam -O ${prefix}.dedup.bam -M ${prefix}.dedup.metrics.txt
## no recalibration unless dealing with human data (low SNP density, many HQ known SNPs)
## call genotype on each sample
/opt/apps/gatk/4.1.9.0/gatk HaplotypeCaller -R $ref -I ${prefix}.dedup.bam --minimum-mapping-quality 30 -ERC GVCF -O ${prefix}.g.vcf.gz
## end merge and call SNPs
