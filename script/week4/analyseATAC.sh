#!/bin/bash
#SBATCH --job-name=makeBigWig      ## Name of the job.
#SBATCH -A ecoevo283         ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-20         ## number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs

module load samtools/1.10
module load ucsc-tools/v393 

inpath="/data/class/ecoevo283/nargesr/ATACseq/"
ref="/data/class/ecoevo283/nargesr/ref/dm6.chrom.sizes"
file=$inpath"prefixes.txt" # prefix file
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`

Nreads=`samtools view -c -F 4 ${inpath}mapped/${prefix}.sort.bam`
Scale=`echo "1.0/($Nreads/1000000)" | bc -l`

samtools view -b ${inpath}mapped/${prefix}.sort.bam | genomeCoverageBed -ibam - -bg -scale $Scale > ${inpath}coverage/${prefix}.coverage

sortBed -i ${inpath}coverage/${prefix}.coverage > ${inpath}coverage/${prefix}.sort.coverage

bedGraphToBigWig ${inpath}coverage/${prefix}.sort.coverage $ref ${inpath}coverage/${prefix}.bw

rm ${inpath}coverage/${prefix}.coverage
