from .cohort import *

def make_cohorts_dict(names,
					file_dirs,
					treats,refs,
					replicates_files,
					sample_groups_files,
					data_dir="../../data/",
					uniprot_file="uniprot-all_20171124.tab.gz"):
	"""
	Create cohort objects via dictionaries of attributes where the key of each dfictionary is the name of the cohort object
	
	Parameters:
	----------


	Returns:
	------

	objs
		Dictionary where keys are each cohort object name and values are the instantiaated cohort object.
	"""
	if len(names)==0:
		return {}
	else:
		objs = {}
		for cohort in names:
			obj = Cohort(cohort=cohort,
								data_dir=data_dir,file_dir=file_dirs[cohort],
								replicates_file=replicates_files[cohort],
								sample_groups_file=sample_groups_files[cohort],
								uniprot_file=uniprot_file)
			obj.set_ref(refs[cohort])
			obj.set_treat(treats[cohort])
			objs[cohort] = obj
		return objs
