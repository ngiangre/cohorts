"""
Testing cohorts package using sample data
"""
from ..cohorts import *

#setting input parameters to make cohorts
top = "~/Github/cohorts/"
file_dirs ={ 'sample' : top+"sample_data/"}
replicates_files = { 'sample' : "df_replicates.tsv"}
sample_groups_files = { 'sample' : "df_sample_groups.tsv"}
name = 'sample'
obj = cohorts.Cohort(cohort = name,
                        file_dir=file_dirs[name],
                        replicates_file=replicates_files[name],
                        sample_groups_file=sample_groups_files[name])

treats = { 'sample' : 'trt' }
refs = {'sample' : 'ref' }
objs = cohorts.cohort.make_cohorts_dict(names=[name],
                        file_dirs=file_dirs,
                        replicates_files=replicates_files,
                        sample_groups_files=sample_groups_files,
                        treats=treats,refs=refs)

objs

print(obj.__doc__)
