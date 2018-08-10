"""
Testing cohorts package using sample data
"""
# Setting working directory
import os
os.chdir('/Users/npg2108/Github/cohorts/src/')

#Importing modules
import sys
sys.path.append("../cohorts")
from utils import *
from cohort import *

#setting input parameters to make cohorts
top = "~/Github/cohorts/"
file_dirs ={ 'sample' : top+"sample_data/"}
replicates_files = { 'sample' : "df_replicates.tsv"}
sample_groups_files = { 'sample' : "df_sample_groups.tsv"}
name = 'sample'
obj = Cohort(cohort = name,
                        file_dir=file_dirs[name],
                        replicates_file=replicates_files[name],
                        sample_groups_file=sample_groups_files[name])


print(obj.__doc__)