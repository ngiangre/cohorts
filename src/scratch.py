"""
For trying out new features
"""

import pkg_resources
pkg_resources.require("cohorts==0.2a0")
import cohorts

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


attribute = 'replicates_hq'
avail_attributes = ['sample_groups', 'replicate_groups', 'raw_replicates', 'raw_samples','replicates_hq','samples_hq']

def get_tidy_data(attribute=""):
    """
    Returns tidy version of available attributes
    
    Parameters
    ----------
    
    attribute: str
        one of ['sample_groups', 'replicate_groups', 'raw_replicates', 'raw_samples','replicates_hq','samples_hq']
    
    Output:
    -------
    
    tidy version of attribute
    
    """
    
    avail_attributes = ['sample_groups', 'replicate_groups', 'raw_replicates', 'raw_samples','replicates_hq','samples_hq']

    if attribute not in avail_attributes:
        raise Exception("Not one of:\n{}".format(avail_attributes))
    elif attribute == 'sample_groups':
        tidy = (getattr(obj,attribute)
            .reset_index()
            .rename(columns={ 'index' : 'Group' })
            .melt(id_vars='Group',var_name='Sample')
            .query('value==1')
            .drop(['value'],axis=1)
            )
    elif attribute == 'replicate_groups':
        tidy = (getattr(obj,attribute)
            .reset_index()
            .rename(columns={ 'index' : 'Group' })
            .melt(id_vars='Group',var_name='Replicate')
            .query('value==1')
            .drop(['value'],axis=1)
            )
    elif attribute == 'raw_replicates':
        tidy = (getattr(obj,attribute)
            .reset_index()
            .rename(columns={ 'index' : 'Protein' })
            .melt(id_vars='Protein',var_name='Replicate',value_name='Value')
            )
    elif attribute == 'replicates_hq' and getattr(obj,attribute) is not None:
        tidy = (getattr(obj,attribute)
            .reset_index()
            .rename(columns={ 'index' : 'Protein' })
            .melt(id_vars='Protein',var_name='Replicate',value_name='Value')
            )
    elif attribute == 'raw_samples':
        tidy = (getattr(obj,attribute)
            .reset_index()
            .rename(columns={ 'index' : 'Protein' })
            .melt(id_vars='Protein',var_name='Sample',value_name='Value')
            )
    elif attribute == 'samples_hq' and getattr(obj,attribute) is not None:
        tidy = (getattr(obj,attribute)
            .reset_index()
            .rename(columns={ 'index' : 'Protein' })
            .melt(id_vars='Protein',var_name='Sample',value_name='Value')
            )
    return tidy

get_tidy_data('sample_groups')