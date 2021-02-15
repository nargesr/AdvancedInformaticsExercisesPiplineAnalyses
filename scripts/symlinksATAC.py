import os
import numpy as np
import pandas as pd

#check directory exist and if not create
directory = "/data/class/ecoevo283/nargesr/ATACseq"
if not os.path.exists(directory):
    os.makedirs(directory)

directory = "/data/class/ecoevo283/nargesr/ATACseq/data"
if not os.path.exists(directory):
    os.makedirs(directory)

cmd = "find /data/class/ecoevo283/public/RAWDATA/ATACseq/*.fq.gz > /data/class/ecoevo283/nargesr/ATACseq/atac_samples.txt"
os.system(cmd)

samples = pd.read_csv("/data/class/ecoevo283/nargesr/ATACseq/atac_samples.txt", sep='\t').values
samples = np.append(samples, np.zeros((samples.shape[0], 1)), axis=1)


for i in range(samples.shape[0]):
        tmp = samples[i, 0].split("4R009_L1_")
        samples[i, 1] = tmp[1]
        cmd = "ln -s " + samples[i, 0] + " /data/class/ecoevo283/nargesr/ATACseq/data/" + samples[i, 1]
        os.system(cmd)
	
pd.DataFrame(samples).to_csv("/data/class/ecoevo283/nargesr/ATACseq/atac_samples.txt", sep='\t')
