"""
Generating sample data to try cohorts package
"""

import numpy as np
import pandas as pd

out_dir = "~/Github/cohorts/sample_data/"

#####helper functions#####
def make_replicates_from_samples(samples,n_replicates=2):

	tmp = [s + "_Rep" for s in samples]
	replicates_bare = np.repeat(tmp,n_replicates)
	order = np.arange(1,n_replicates + 1,1)
	rep_nums = np.tile(order,n_samples)
	replicates = []
	for i, num in enumerate(rep_nums):
	    replicates.append(replicates_bare[i] + str(num))

	return replicates


def make_df_replicates(markers,replicates,loc1=2,loc2=20):

	mat_num = len(markers) * len(replicates)
	marker_values_ref = np.random.exponential(loc1,mat_num//2)
	marker_values_trt = np.random.exponential(loc2,mat_num//2)
	df_replicates_ref = pd.DataFrame(marker_values_ref.reshape(len(markers),len(replicates)//2),
	                                index = markers, 
	                                columns = replicates[:len(replicates)//2])
	df_replicates_trt = pd.DataFrame(marker_values_trt.reshape(len(markers),len(replicates)//2),
	                                index = markers, 
	                                columns = replicates[len(replicates)//2:len(replicates)])
	df_replicates = df_replicates_ref.join(df_replicates_trt)

	return df_replicates

def make_df_samples(markers,samples,loc1=2,loc2=20):

	mat_num = len(markers) * len(samples)
	marker_values_ref = np.random.exponential(loc1,mat_num//2)
	marker_values_trt = np.random.exponential(loc2,mat_num//2)
	df_samples_ref = pd.DataFrame(marker_values_ref.reshape(len(markers),len(samples)//2),
	                                index = markers, 
	                                columns = samples[:len(samples)//2])
	df_samples_trt = pd.DataFrame(marker_values_trt.reshape(len(markers),len(samples)//2),
	                                index = markers, 
	                                columns = samples[len(samples)//2:len(samples)])
	df_samples = df_samples_ref.join(df_samples_trt)

	return df_samples

def make_df_sample_groups(samples):

	n_samples = len(samples)
	trt_status = np.repeat([1,0],np.floor(n_samples/2).astype(int))
	ref_status = np.repeat([0,1],np.ceil(n_samples/2).astype(int))
	n_cov = 2
	tmp = np.repeat(1,n_cov)
	cov_status = np.pad(tmp,(0,n_samples-n_cov),mode="constant")
	df_sample_groups = pd.DataFrame( [ trt_status, ref_status,cov_status ],
	                columns = samples,index = ["trt","ref","cov"])

	return df_sample_groups

##########	

#Set markers - here, I picked the marker to be proteins, so I just picked random Uniprot identifiers
proteins = [ 'E7EX29', 'P62191', 'Q99460', 'P52209', 
'P08253', 'P27338', 'P15144', 'P05067', 'Q9UJX5', 
'Q8N302', 'P03950', 'Q9Y5C1']

#####SAMPLE DATA#####

#Make Samples
n_samples = 10
samples = ["S"+ str(i) for i in np.arange(1,n_samples + 1,1)]

#Make Replicates
replicates = make_replicates_from_samples(samples)

#Make df_replicates dataframes
df_replicates = make_df_replicates(markers=proteins,replicates=replicates,loc1=2,loc2=20)
df_replicates.to_csv(out_dir+"df_replicates.tsv",sep="\t")

#Make df_samples dataframe
df_samples = make_df_samples(markers=proteins,samples=samples,loc1=2,loc2=20)
df_samples.to_csv(out_dir+"df_samples.tsv",sep="\t")

#Make df_sample_groups dataframe
df_sample_groups = make_df_sample_groups(samples)
df_sample_groups.to_csv(out_dir+"df_sample_groups.tsv",sep="\t")

##########


#####PATIENT COHORT 1#####

out_dir = "~/Github/cohorts/sample_data/patient_cohort_1/"

#Make Samples
n_samples = 10
samples = ["P1-"+ str(i) for i in np.arange(1,n_samples + 1,1)]

#Make Replicates
replicates = make_replicates_from_samples(samples)

#Make df_replicates dataframes
df_replicates = make_df_replicates(markers=proteins,replicates=replicates,loc1=2,loc2=20)
df_replicates.to_csv(out_dir+"df_replicates.tsv",sep="\t")

#Make df_samples dataframe
df_samples = make_df_samples(markers=proteins,samples=samples,loc1=2,loc2=20)
df_samples.to_csv(out_dir+"df_samples.tsv",sep="\t")

#Make df_sample_groups dataframe
df_sample_groups = make_df_sample_groups(samples)
df_sample_groups.to_csv(out_dir+"df_sample_groups.tsv",sep="\t")

#####PATIENT COHORT 1#####

out_dir = "~/Github/cohorts/sample_data/patient_cohort_2/"

#Make Samples
n_samples = 10
samples = ["P2-"+ str(i) for i in np.arange(1,n_samples + 1,1)]

#Make Replicates
replicates = make_replicates_from_samples(samples)

#Make df_replicates dataframes
df_replicates = make_df_replicates(markers=proteins,replicates=replicates,loc1=3,loc2=21)
df_replicates.to_csv(out_dir+"df_replicates.tsv",sep="\t")

#Make df_samples dataframe
df_samples = make_df_samples(markers=proteins,samples=samples,loc1=3,loc2=21)
df_samples.to_csv(out_dir+"df_samples.tsv",sep="\t")

#Make df_sample_groups dataframe
df_sample_groups = make_df_sample_groups(samples)
df_sample_groups.to_csv(out_dir+"df_sample_groups.tsv",sep="\t")

#####PATIENT COHORT 1#####

out_dir = "~/Github/cohorts/sample_data/patient_cohort_3/"

#Make Samples
n_samples = 10
samples = ["P3-"+ str(i) for i in np.arange(1,n_samples + 1,1)]

#Make Replicates
replicates = make_replicates_from_samples(samples)

#Make df_replicates dataframes
df_replicates = make_df_replicates(markers=proteins,replicates=replicates,loc1=1,loc2=19)
df_replicates.to_csv(out_dir+"df_replicates.tsv",sep="\t")

#Make df_samples dataframe
df_samples = make_df_samples(markers=proteins,samples=samples,loc1=1,loc2=19)
df_samples.to_csv(out_dir+"df_samples.tsv",sep="\t")

#Make df_sample_groups dataframe
df_sample_groups = make_df_sample_groups(samples)
df_sample_groups.to_csv(out_dir+"df_sample_groups.tsv",sep="\t")
