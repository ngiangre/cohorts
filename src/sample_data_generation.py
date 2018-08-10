"""
Generating sample data to try cohorts package
"""

import numpy as np
import pandas as pd

out_dir = "~/Github/cohorts/sample_data/"

#Make Samples
n_samples = 10
samples = ["S"+ str(i) for i in np.arange(1,n_samples + 1,1)]

#Make Replicates
tmp = [s + "_Rep" for s in samples]
n_replicates = 2
replicates_bare = np.repeat(tmp,n_replicates)
order = np.arange(1,n_replicates + 1,1)
rep_nums = np.tile(order,n_samples)
replicates = []
for i, num in enumerate(rep_nums):
    replicates.append(replicates_bare[i] + str(num))

#Set Proteins - just picked random Uniprot identifiers
proteins = [ 'E7EX29', 'P62191', 'Q99460', 'P52209', 
'P08253', 'P27338', 'P15144', 'P05067', 'Q9UJX5', 
'Q8N302', 'P03950', 'Q9Y5C1']

#Make df_replicates dataframe
mat_num = len(proteins) * len(replicates)
protein_values = np.random.exponential(2,mat_num)
df_replicates = pd.DataFrame(protein_values.reshape(len(proteins),len(replicates)),
                                index = proteins, columns = replicates)
df_replicates.to_csv(out_dir+"df_replicates.tsv",sep="\t")

#Make df_sample_groups dataframe
trt_status = np.repeat([1,0],np.floor(n_samples/2).astype(int))
ref_status = np.repeat([0,1],np.ceil(n_samples/2).astype(int))
n_cov = 2
tmp = np.repeat(1,n_cov)
cov_status = np.pad(tmp,(0,n_samples-n_cov),mode="constant")
df_sample_groups = pd.DataFrame( [ trt_status, ref_status,cov_status ],
                columns = samples,index = ["trt","ref","cov"])
df_sample_groups.to_csv(out_dir+"df_sample_groups.tsv",sep="\t")


