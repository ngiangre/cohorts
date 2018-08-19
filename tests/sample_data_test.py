from unittest import TestCase

import cohorts
from os import getcwd

class TestCohort(TestCase):
	def test_cohort_name_is_string(self):

		top = getcwd() + "/"
		print(top)
		file_dirs ={ 'sample' : top+"sample_data/"}
		replicates_files = { 'sample' : "df_replicates.tsv"}
		sample_groups_files = { 'sample' : "df_sample_groups.tsv"}
		name = 'sample'
		obj = cohorts.Cohort(cohort = name,
							file_dir=file_dirs[name],
							replicates_file=replicates_files[name],
							sample_groups_file=sample_groups_files[name])
		self.assertTrue(isinstance(obj.cohort, str))