import os
import numpy as np
import pandas as pd

#check directory exist and if not create
directory = "/data/class/ecoevo283/nargesr/DNAseq"
if not os.path.exists(directory):
    os.makedirs(directory)

cmd = "find /data/class/ecoevo283/public/RAWDATA/DNAseq/*.fq.gz > /data/class/ecoevo283/nargesr/DNAseq/dna_samples.txt"
os.system(cmd)

samples = pd.read_csv("/data/class/ecoevo283/nargesr/DNAseq/dna_samples.txt", sep='\t').values
samples = np.append(samples, np.zeros((samples.shape[0], 1)), axis=1)

for i in range(samples.shape[0]):
	tmp = samples[i, 0].split("q/")
	samples[i, 0] = tmp[0]
	samples[i, 1] = tmp[1]
	cmd = "ln -s " + samples[i, 1] + " /data/class/ecoevo283/nargesr/DNAseq/data/" + samples[i, 1]
	os.system(cmd)

pd.DataFrame(samples).to_csv("/data/class/ecoevo283/nargesr/DNAseq/dna_samples.txt", sep='\t')
