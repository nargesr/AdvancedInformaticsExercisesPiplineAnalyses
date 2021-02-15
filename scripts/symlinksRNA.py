import os
import numpy as np
import pandas as pd

#check directory exist and if not create
directory = "/data/class/ecoevo283/nargesr/RNAseq"
if not os.path.exists(directory):
    os.makedirs(directory)

directory = "/data/class/ecoevo283/nargesr/RNAseq/data"
if not os.path.exists(directory):
    os.makedirs(directory)


samples = pd.read_csv("/data/class/ecoevo283/public/RAWDATA/RNAseq/RNAseq384_SampleCoding.txt", sep='\t')

samples['sampleName'] = samples.apply(lambda x: "x" + x.FullSampleName, axis=1)
samples['plex'] = samples.apply(lambda x: x.Multiplexi5index.split("_")[0], axis=1)

pd.DataFrame(samples).to_csv("/data/class/ecoevo283/nargesr/RNAseq/rna_samples.txt", sep='\t')


for i in range(samples.shape[0]):
	path = "/data/class/ecoevo283/public/RAWDATA/RNAseq/RNAseq384plex_flowcell01/Project_" + samples['plex'][i] + "/Sample_" + str(samples['SampleNumber'][i])
	cmd = "ln -s " + path + "/*_R1_001.fastq.gz" + " /data/class/ecoevo283/nargesr/RNAseq/data/" + samples['sampleName'][i] + "_R1.fq.gz"
	os.system(cmd)
	cmd = "ln -s " + path + "/*_R2_001.fastq.gz" + " /data/class/ecoevo283/nargesr/RNAseq/data/" + samples['sampleName'][i] + "_R2.fq.gz"
	os.system(cmd)

