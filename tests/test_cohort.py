'''
import unittest
from cohorts import Cohort

class TestCohortMethods(object):


	def __init__(self):

		top = "./"
		data_dirs = { 'my_cohort' : top+"data/" }
		file_dirs = { 'my_cohort' : top+"sample_data/" }
		replicates_files = { 'my_cohort' : "df_replicates.tsv" }
		#sample_files = { 'my_cohort' : "df_samples.tsv" }
		sample_groups_files = { 'my_cohort' : "df_sample_groups.tsv" }
		cohort = { 'my_cohort' : 'my_cohort' }
		references = { 'my_cohort' : 'ref' }
		treatments = { 'my_cohort' : 'trt' }
		cohort_name = 'my_cohort'
		c = Cohort(cohort=cohort[cohort_name],marker_type="protein",
							data_dir=data_dirs[cohort_name],file_dir=file_dirs[cohort_name],
							replicates_file=replicates_files[cohort_name],
							#sample_files=sample_files[cohort],
							sample_groups_file=sample_groups_files[cohort_name],
							reference=references[cohort_name],
							treatment=treatments[cohort_name])

		self.cohort = c

	def test_cohort_name(self):
		return self.cohort.cohort=='my_cohort'

if __name__ == '__main__':
	unittest.main()
'''
def given():
	return True