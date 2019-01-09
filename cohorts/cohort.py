from .utils import *

import os
import numpy as np
import pandas as pd
import scipy.stats as sc
from sklearn.preprocessing import StandardScaler
import itertools as it

class Cohort(object):
	"""
	Patient cohort object for patient proteomics data.

	Dataframes, variables, and functions facilitating the processing, analysis, and integration of the cohort data.

	Parameters
	----------
	cohort: str
		name of the patient cohort

	file_dir: str
		directory where replicate dataframe and 
		sample group membership dataframe file names are located

	replicates_file: str
		name of the replicates dataframe file

		A markers x B replicates
		comma (*.csv) or tab (*.tsv) delimited

		replicate = "SampleName" + "_Rep[0-9]"

	sample_groups_file: str
		name of the sample group file
		comma (*.csv) or tab (*.tsv) delimited

		N groups x M samples

		sample = "SampleName"

	data_dir: str
		directory where extra data files are located
		
	uniprot_file: str
		name of the uniprot database flat file located in data_dir

	Examples
	--------
	>>>c = cohorts.Cohort(cohort='cohort_name',
		file_dir="path/to/files/",
		replicates_file="file_name",
		samples_file="file_name",
		sample_groups_file="sample_groups_file_name",
		data_dir="path/to/data/dir/",
		marker_type="protein",
		reference="ref",treatment="trt"
		)
	>>>c.set_replicates_hq()
	>>>c.set_trans_replicates_hq()
	>>>c.set_samples_hq()
	>>>c.set_trans_samples_hq()

	"""

	tests = [ 
					( "t-test",sc.ttest_ind ),
					("Wilcoxon_RankSum_test",sc.ranksums)
				]


	def __init__(self,cohort='cumc',
		replicates_file=None,
		samples_file=None,
		sample_groups_file=None,
		file_dir="",
		data_dir="../../data/",
		marker_type="protein",
		reference=None,
		treatment=None):

		# OBJECT CORE ATTRIBUTES
		self.cwd = os.getcwd()
		self.data_dir = data_dir
		self.file_dir = file_dir
		self.cohort = cohort
		self.marker_type = marker_type
		
		self.reference = reference
		self.treat = treatment
		self.data = {}
		self.params = {}

		# REPLICATE SPECIFIC ATTRIBUTES
		self.replicates_file = replicates_file
		self.raw_replicates = None
		self.replicates_hq = None
		self.trans_replicates_hq = None
		self.sample_replicate_dictionary = None
		self.replicate_groups = None
		self.replicates = None
		
		# SAMPLE SPECIFIC ATTRIBUTES
		self.samples_file = samples_file
		self.sample_groups_file = sample_groups_file
		self.raw_samples = None
		self.samples_hq = None
		self.trans_samples_hq = None
		self.samples = None
		self.sample_groups = None
		self.groups = None
		self.markers = None
		
		# SET ATTRIBUTE FUNCTIONS
		if replicates_file is not None:
			self.set_replicates_file()
			self.set_raw_replicates()
			self.set_replicates()
			self.set_sample_replicate_dictionary()
			self.set_raw_samples(got_replicates=True)

			self.set_sample_groups_file()
			self.set_sample_groups()
			self.set_groups()
			self.set_samples()
			self.set_markers()
				
			self.set_replicate_groups()

		
	#SET FUNCTIONS
	
	def set_replicates_file(self):
		"""
		Setting the replicates file string.
		Combination of the file directory string and the replicates file string

		A csv file is needed!
		Parameters
		----------
		None

		"""

		self.replicates_file = str(self.file_dir) + str(self.replicates_file)

	def set_sample_groups_file(self):
		"""
		Setting the replicate_groups string
		Combination of the path string and the replicates file string

		A csv file is needed!
		
		Parameters
		----------
		None

		"""
		self.sample_groups_file = str(self.file_dir) + str(self.sample_groups_file)

	def set_raw_replicates(self):
		"""
		Setting raw replicate dataframe from files given to object

		Parameters
		----------
		None

		"""		

		#checking if tab separated (tsv) or comma separated (csv) format
		splt = self.replicates_file.split(".")
		format = splt[len(splt)-1][0]
		if format is 'c':
			self.raw_replicates = pd.read_csv(self.replicates_file,delimiter=",",index_col=0).sort_index(axis=0).sort_index(axis=1)
		else:
			self.raw_replicates = pd.read_csv(self.replicates_file,delimiter="\t",index_col=0).sort_index(axis=0).sort_index(axis=1)
		
		self.raw_replicates.index.name = self.marker_type
		self.raw_replicates.columns.name = "replicate"

	def set_replicates(self):
		"""
		Setting replicate names 

		Depends on : raw_replicates

		Parameters
		----------
		None

		"""

		self.replicates = self.raw_replicates.columns.tolist()

	def set_markers(self):
		"""
		Setting marker names 

		Depends on : raw_replicates

		Parameters
		----------
		None

		"""

		self.markers = self.raw_samples.index.tolist()

	def set_sample_replicate_dictionary(self):
		"""
		Setting dictionary of samples (keys) to replicates (values)

		Depends on : replicates

		Parameters
		----------
		None

		"""

		#get replicates
		replicates = self.replicates

		#strip replicate indicator to get array of sample names for replicates
		duplicate_samples = [y[0] for y in [x for x in np.chararray.split(replicates,"_")]]

		#set previous to series for easy unique sample name retrieval
		samples = pd.Series(duplicate_samples).unique()

		#initialize dictionary
		dictionary = {}

		#loop through unique samples
		for key in samples:
			#get indice of replicates for sample
			replicate_inds = pd.Series(duplicate_samples).isin([key])

			#get replicates for sample
			values = np.array(replicates)[replicate_inds.tolist()]

			#set key:value pair (sample : replicates) in dictionary
			dictionary[key] = values.tolist()

		self.sample_replicate_dictionary = dictionary

	def set_raw_samples(self,agg='mean',got_replicates=True):
		"""
		Setting raw sample dataframe 

		Take statistic of replicates with particular sample membership

		Depends on : raw_replicates and make_df_samples()

		Parameters
		----------

		agg: {'mean','median', 'variance'}
			statistic to apply to replicate values for sample membership

		"""

		if got_replicates:
			df = self.raw_replicates.fillna(0)
			
			self.raw_samples = self.make_df_samples(df,agg=agg)
		else:
			df = self.raw_replicates.fillna(0)
			
			self.raw_samples = self.make_df_samples(df,agg=agg)

	def set_replicates_hq(self,
		uniprot_annot=False,
		annot_status='reviewed',
		quant_least_reps_per_samps=False,
		n_reps=1,
		all_reps_quant=False,
		threshold=88,
		intrasamp_var=False,
		higher=False):
		"""
		Setting markers in raw dataframe to be reviewed by Uniprot

		Depends on : raw_replicates, markers_in_n_replicates_per_samples, markers_quant_all_reps, proteins_with_uniprot_annot, markers_by_intrasample_variability

		Parameters
		----------
		
		uniprot_annot {True,False}
			Subset proteins by annot=ation status in uniprot database
		annot_status {'reviewed','unreviewed'}
			Annotation status in uniprot
		quant_least_reps_per_samps {True,False}
			Subset proteins by the number of replicates quantified per sample
		n_reps
			Number of replicates quantified per sample. See above
		all_reps_quant {True,False}
			Subset markers by only those quantified in every replicate per sample
		intrasamp_var {True,False}
			Subset markers by the amount of replicate variance per sample
		threshold
			Quantile to subset variance from. See above
		higher {True,False}
			Indicates whether to subset by proteins with a high annotation score from Uniprot

		"""

		#get raw data
		df = self.raw_replicates.copy().fillna(0)

		#subset to have markers found in atleast sufficient_reps 
		#per sample for all samples
		if quant_least_reps_per_samps:
			markers = self.markers_in_n_replicates_per_samples(df,n_reps=n_reps)
			df = df.loc[markers]

		#subset to have proteins found in atleast n_reps 
		#per sample for all samples
		if all_reps_quant:
			markers = self.markers_quant_all_reps(df)
			df = df.loc[markers]

		# subset to have markers with certain intrasample variability
		if intrasamp_var:

			markers = self.markers_by_intrasample_variability(df,
							higher=higher,
							threshold=threshold)
			df = df.loc[markers]

		self.replicates_hq = df

	def set_trans_replicates_hq(self,
		trans='None',
		add_small=False,
		stat_thresh=False,
		threshold=95,
		higher=False,
		statistic='variance'):
		"""
		Transforming hq values with different functions, by default log1p.
		
		Depends on : replicates_hq, markers_by_statistic_threshold

		Parameters
		----------

		trans: numerical transformation. Default: scikitlearn's StandardScaler
			Indicates how to transform the raw marker values
		add_small
			Add small value if any dataframe value is 0 - when log transforming
		stat_threshold
			Subsetting dataframe by a statistic threshold
		threshold [0,100]
			Quantile for statistic threshold
		higher
			Subset above or below statistic threshold
		statistic {'mean','median','variance'}
			Statistic to be applied

		"""

		#set hq dataframe; make sure there's no NaN values
		df = self.replicates_hq

		#instantiate/declare/set up Standard scaler model from
		#scikitlearn on data
		scaler = StandardScaler()
		scaler.fit(df)

		#list of name-function pairs
		func = [('log1p', np.log1p),
				('log2', np.log2),
				('sklearn', scaler.transform),
				('rank_normalize', rank_normalize),
				('quantile_normalize', quantileNormalize),
				('None', pd.DataFrame.copy)]

		#get appropriate list index for function from trans string
		m = [x for x in range(len(func)) if trans in func[x][0]][0]

		#add small epsilon value to all protein values
		if add_small:
			eps = 0.0001
			df = df.applymap(lambda x : x + eps)

		#apply function from appropriate function list index
		data = pd.DataFrame(func[m][1](df),
								index=df.index,
								columns=df.columns
							)

		#subset to have proteins that meet a certain statistical threshold
		if stat_thresh:
			markers = self.markers_by_statistic_threshold(data,
							statistic=statistic,
							threshold=threshold,
							higher=higher)
			data = data.loc[markers]

		self.trans_replicates_hq = data

	def set_samples_hq(self,
		agg='mean',
		uniprot_annot=False,
		annot_status='reviewed'):
		"""
		Setting sample dataframe from processed replicate dataframe

		Depends on : trans_replicates_hq, make_df_samples
		Parameters
		----------
		
		agg: {'mean','median', 'variance'}
			statistic to apply to replicate values for sample membership
		uniprot_annot {True,False}
			Subset proteins by annot=ation status in uniprot database
		annot_status {'reviewed','unreviewed'}
			Annotation status in uniprot
		"""

		#get processed replicates
		df = self.trans_replicates_hq

		self.samples_hq = self.make_df_samples( df , agg = agg )

	def set_trans_samples_hq(self,
		trans='None',
		add_small=False,
		stat_thresh=False,
		threshold=95,
		higher=False,
		statistic='variance'):
		"""
		Transforming hq values with different functions, by default log1p.
		
		Depends on : samples_hq, markers_by_statistic_threshold

		Parameters
		----------

		trans: numerical transformation. Default: scikitlearn's StandardScaler
			Indicates how to transform the raw marker values
		add_small
			Add small value if any dataframe value is 0 - when log transforming
		stat_threshold
			Subsetting dataframe by a statistic threshold
		threshold [0,100]
			Quantile for statistic threshold
		higher
			Subset above or below statistic threshold
		statistic {'mean','median','variance'}
			Statistic to be applied

		"""

		#set hq dataframe; make sure there's no NaN values
		data = self.samples_hq

		#instantiate/declare/set up Standard scaler model from
		#scikitlearn on data
		scaler = StandardScaler()
		scaler.fit(data)

		#only do log (other than log1p) transformations if values are
		#all nonzero
		if np.sum((data==0).values)>0:
			add_small = True
			print('There are zero values in the sample dataframe, a small number epsilon must be added to the marker values when there are zero values. Adding small epsilon...')

		#list of name-function pairs
		func = [('log1p', np.log1p),
				('log2', np.log2),
				('sklearn', scaler.transform),
				('rank_normalize', rank_normalize),
				('quantile_normalize', quantileNormalize),
				('None', pd.DataFrame.copy)]

		#get appropriate list index for function from trans string
		m = [x for x in range(len(func)) if trans in func[x][0]][0]

		#add small epsilon value to all protein values
		if add_small:
			eps = 0.0001
			data = data.applymap(lambda x : x + eps)

		#apply function from appropriate function list index
		data = pd.DataFrame(func[m][1](data),
								index=data.index,
								columns=data.columns
							)

		#subset to have proteins that meet a certain statistical threshold
		if stat_thresh:
			markers = self.markers_by_statistic_threshold(data,
							statistic=statistic,
							threshold=threshold,
							higher=higher)
			data = data.loc[markers]

		self.trans_samples_hq = data

	def set_samples(self):
		"""
		Setting patient names from raw dataframe

		Depends on : sample_replicate_dictionary

		Parameters
		----------
		None

		"""
		self.samples = [x for x in self.sample_replicate_dictionary.keys()]

	def set_sample_groups(self,file=None):
		"""
		Setting group membership of samples dataframe where 1 indicates membership and 0 no membership. 

		Parameters
		----------
		file
			file to read sample groups data. Defaults to parameters from cohorts instantiation
		"""		

		#checking if tab separated (tsv) or comma separated (csv) format
		splt = self.sample_groups_file.split(".")
		format = splt[len(splt)-1][0]
		if format is 'c':
			self.sample_groups = pd.read_csv(self.sample_groups_file,delimiter=",",index_col=0).sort_index(axis=0).sort_index(axis=1)
		else:
			self.sample_groups = pd.read_csv(self.sample_groups_file,delimiter="\t",index_col=0).sort_index(axis=0).sort_index(axis=1)
			
		self.sample_groups.columns.name = "sample"
		self.sample_groups.index.name = "group"

	def set_df_replicate_groups(self):
		"""
		Aggregating df_replicate_groups to derive sample group 
		membership

		Parameters
		----------
		None

		"""

		mats = {}
		for samp in self.samples:
			single = self.sample_groups.loc[:,samp]
			reps = self.sample_replicate_dictionary[samp]
			for r in reps:
				mats[r] = single

		df = pd.DataFrame.from_dict(mats)

		df.columns.name = "replicate"
		df.index.name = "group"

		return df

	def set_replicate_groups(self,file=None):
		"""
		Setting group membership of replicate dataframe where 1 indicates membership and 0 no membership. 
		
		Parameters
		----------
		file
			file to read replicate groups data. Defaults to parameters from cohorts instantiation

		"""
			
		self.replicate_groups = self.set_df_replicate_groups()

	def set_groups(self):
		"""
		Setting groups from sample groups dataframe
		
		Parameters
		----------
		None
		
		"""

		self.groups = self.sample_groups.index.values
 
	def set_ref(self,ref='NL'):
		"""
		Setting reference group 
		
		Parameters
		----------
		
		ref: string
			Reference group name
		
		
		"""

		#set new reference variable if it's in the groups already. If not, just ignore it. 
		if np.any(self.groups==ref):

			self.reference = ref

	def set_treat(self,trt='PGD'):
		"""
		Setting treatment group
		
		Parameters
		----------
		
		trt: string
			Treatment group name
		"""

		#set new reference variable if it's in the groups already. If not, just ignore it. 
		if np.any(self.groups==trt):

			self.treat = trt

	def get_tidy_data(self,attribute=""):
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
		
		attributes = [x for x in self.__dict__.keys()]
		attribute_types = [type(getattr(self,x)) for x in self.__dict__.keys()]
		inds = [x==pd.core.frame.DataFrame for x in attribute_types]
		avail_attributes = np.array(attributes)[np.array(inds)].tolist()

		if attribute not in avail_attributes:
			raise Exception("Not one of:\n{}".format(avail_attributes))

		df = getattr(self,attribute).copy()

		tidy = df.reset_index().melt(id_vars=df.index.name)

		return tidy


	#ANALYSIS FUNCTIONS

	def manual_feature_extraction(self,df):
		"""
		Query protein quantification between reference and other groups

		Parameters
		----------
		df
			dataframe to manually extract features
		"""
		
		#load raw and sample group dataframes
		df_samples = df
		df_sample_groups = self.sample_groups.copy()

		#make rownames-presence/absence/mixed conditions of markers amongst samples
		val_grps = ('allq', 'allnotq', 'mixed')
		tmp = list(it.product(val_grps,val_grps))
		rownames = []
		for i in tmp:
			rownames.append("_".join(i))

		#set column names-reference/indicator scenarios
		x = self.groups != self.reference
		t = [ self.reference + '/' + x for x in self.groups[x] ]
		colnames = tuple(t)
		

		#populate length and marker array dataframes
		df_len = pd.DataFrame(index=rownames,columns=colnames)
		df_arr = pd.DataFrame(index=rownames,columns=colnames)

		#populate each element with an empty array or list
		for i, j in it.product(range(df_len.shape[0]),range(df_len.shape[1])):
			df_len.iloc[i][j] = ()
		for i, j in it.product(range(df_arr.shape[0]),range(df_arr.shape[1])):
			df_arr.iloc[i][j] = []

		#make list of functions for assessing presence/absence/mixed protein status
		singlefuncs = [ ( "allq", allq ) , 
						("allnotq", allnotq ) ,
						("mixed", mixed ) 
					]

		#loop through indicators to assess presense/absence/mixed proteins for each reference/indicator scenario and populate length/protein array dataframes
		for i in self.groups[x]:

			#set reference/indicator scenario
			gr1 = self.reference
			gr2 = i
			colname = gr1+"/"+gr2
			print(colname)
			
			#make arrays of proteins for reference samples for each assessment condition 
			gr1_arr = []
			tmp = df_samples.T
			arr1 = df_sample_groups.loc[gr1] == 1
			X = tmp[arr1.values].transpose()
			for l in range(len(val_grps)):
				tmp = singlefuncs[l][1](X)#slow
				gr1_arr.append(tmp)
		
			#make arrays of markers for reference samples for each assessment condition 
			gr2_arr = []
			tmp = df_samples.T
			arr2 = df_sample_groups.loc[gr2] == 1
			X = tmp[arr2.values].transpose()
			for k in range(len(val_grps)):
				tmp = singlefuncs[k][1](X)
				gr2_arr.append(tmp)
			
			#for each combination of gr1_arr and gr2_arr, do intersection and append to corresponding row in column of data frame. This gives proteins in each condition for each reference/indicator scenario
			for m,n in it.product( range( len(gr1_arr) ) , range( len(gr2_arr) ) ):
				rowname = singlefuncs[m][0]+"_"+singlefuncs[n][0]
				df_len.loc[ rowname , colname ] = len ( np.intersect1d( gr1_arr[m], gr2_arr[n] ) ) 
				df_arr.loc[ rowname , colname ] = np.intersect1d( gr1_arr[m], gr2_arr[n] ) 

		#dictionary of length and protein array dataframes
		dictionary = { 'df_len' : df_len, 'df_arr' : df_arr} 

		helper_dictionary = self.make_marker_substraction(dictionary)

		self.data['mfe'] = { 'main' : dictionary, 'helper' : helper_dictionary } 

	def make_marker_substraction(self,dictionary=None):
		"""
		Helper function only for manual_feature_extraction method to do set operations on marker results

		Right now supports difference of markers

		Parameters
		----------
		dictionary
			dictionary to fill in-fed from manual_feature_extraction

		Output
		------
		Dictionary of marker array length band array from difference of manually extracted features between ref and treat groups
		"""

		#declare dictionary from manual feature extraction
		mfe = dictionary

		#set comparisons made in mfe
		comps_grps = mfe['df_len'].columns.tolist()

		#set marker feature comparisons from manual feature extraction
		comps_markers = mfe['df_len'].index.tolist() 

		#make combination of comparisons-resulting dataframe rows and columns
		combs = list(
			it.product(
				comps_grps,
				comps_grps
			)
		)

		#make dictionary of presence,absence, mixed. Values will be NxN dataframe labeled with combs above
		markers_dict = dict.fromkeys(
			comps_markers
		)

		#set dataframe to extract markers from
		df = mfe['df_arr']

		#loop through dictionary keys which are the marker feature comparisons
		for i in markers_dict.keys():

			#make empty NxN dataframe to store protein lengths and arrays
			empty_len = pd.DataFrame(
					index = comps_grps,
			 		columns = comps_grps
					)

			#make copied empty dataframe to store protein arrays
			empty_diff = empty_len.copy()

			#loopp through dataframe index, columns and store length of difference and proteins from difference of sets
			for j in comps_grps:
				for k in comps_grps:
					diff_arr = set(df.loc[i,j]) - set(df.loc[i,k])
					length = len(diff_arr)
					empty_len.at[j,k] = length
					empty_diff.at[j,k] = diff_arr

			#store dictionary of dataframes as value in dictionary
			markers_dict[i] = { 'df_len' : empty_len, 'df_arr' : empty_diff }

		return markers_dict

	def hypothesis_testing(self,df,df_groups,tests=None):
		"""
		Hypothesis testing of reference sample proteins versus treatment sample markers.
		
		Parameters
		----------
		df
			dataframe of values to test
		df_groups
			group membership of samples/replicates in df to test

		Outputs
		-------
		df
			hypothesis test results
		"""

		#get subset dataframes by reference and treatment groups
		df1 = self.get_sub_df(df,df_groups,self.reference)
		df2 = self.get_sub_df(df,df_groups,self.treat)

		#List of hypothesis test names and functions
		if tests is None:
			tests = self.tests
		else:
			print('A list of statistical tests are required. For example: [( "t-test",scipy.stats.ttest_ind )]. More than one can be given.')
			tests = tests

		#set dataframe to fill
		df = pd.DataFrame(columns=[self.marker_type,'test','pvalue','statistic'])

		#for each hypothesis test name and function
		for hyp in tests:

			#set test name
			test = hyp[0]

			#for each protein
			for marker in df1.index:

				#get protein location in dataframes
				firstloc = df1.loc[marker]
				secondloc = df2.loc[marker]

				#make sure the arrays are filled with floats
				a = firstloc.values.astype(float)
				b = secondloc.values.astype(float)

				#set statistic from test
				stat = hyp[1](a,b)[0]

				#set pvalue from test
				pval = hyp[1](a,b)[1]

				#put data in pandas series
				row = pd.Series([self.marker_type,test,pval,stat],index=df.columns)

				#add row to dataframe
				df = df.append(row,ignore_index=True)

		#set correct object membership of dataframe columns
		df = df.astype(dtype= {
								self.marker_type:"str",
								"test":"str",
								"pvalue":"float64",
								"statistic":'float64'
							}
						)

		return df

	def get_sub_df(self,df,df_groups,grp):
		"""
		Helper function in class to get patient-subset dataframe.

		Subset larger dataframe by reference and treatment groups

		Used in hypothesis_testing()

		Parameters
		----------
		df
			'omics data dataframe given to hypothesis_testing
		df_groups
			groups dataframe given to hypothesis testing
		grp str
			group name for subsetting

		Output
		------
		subsetted dataset
		"""

		#get indices of replicates/samples in groups
		inds = np.where(df_groups.loc[grp] == 1)

		#get replicates/samples
		samples = df_groups.columns[inds]
		
		#return
		
		return df.loc[:,samples]

	#MARKER SUBSET FUNCTIONS

	def markers_in_n_replicates_per_samples(self,df,n_reps=1,value=0):
		"""
		Subsetting dataset by markers that are quantified in more than n replicates per sample for all samples

		Parameters
		----------
		df
			dataframe to subset
		n_reps
			Least number of replicates quantified per sample

		Output
		------
		subsetted markers
		"""
		samp_rep_dict = self.sample_replicate_dictionary
		
		atleast1_reps = []
		for samp in samp_rep_dict.keys():
			reps = samp_rep_dict[samp]
			sub = df.loc[:,reps]
			mask = sub.apply(lambda x : x > value,axis=1)
			#markers quantified in more than 1 replicates of a sample
			markers_bool = mask.sum(axis=1) > n_reps
			[atleast1_reps.append(x) for x in markers_bool.index[markers_bool].tolist()]
		
		atleast1_rep_per_samp_markers = np.unique(atleast1_reps)

		return atleast1_rep_per_samp_markers

	def markers_quant_all_reps(self,df):
		"""
		Helper function returning all markers that 
		are quantified in every replicate.

		Parameters
		----------
		df
			dataframe to subset

		Output
		------
		subsetted markers
		"""

		#determine how many reps show quantification of the marker
		num_reps_quant = df.apply(lambda x : sum(x > 0),axis=1)

		#is the num reps all of them?
		all_reps_quant = num_reps_quant.apply(
			lambda x : x == df.shape[1])

		#get protein names
		allquant_markers = df.index[all_reps_quant]

		return allquant_markers

	def markers_by_statistic_threshold(self,
		df,
		statistic='mean',
		higher=False,
		threshold=80):
		"""
		Markers that meet a certain statistical threshold
		e.g. markers that have low variance across the dataset
		
		outputs markers above/below that threshold

		e.g. threshold=25 gives markers that are below the 25th quantile
		of the statistical function. Threshold=80 gives proteins that are 
		below the 80th quantile of the statistical function
		
		Parameters
		----------
		df
			dataframe to subset
		statistic: {'mean','median', 'variance'}
			statistic to apply across values
		higher {True,False}
			Indicates whether to subset by markers higher (True) than or less than/equal to (False) the threshold score 
		threshold [0,100]
			Quantile for statistic threshold

		Output
		------
		subsetted markers
		
		"""
		#possible statistical functions to threshold by
		funcs = {
			'mean' : np.mean,
			'median' : np.median,
			'variance' : np.var,
			'std' : np.std
		}
		
		#collect proteins in dataset
		markers = df.index
		
		#create series by applying statistic
		series = df.apply(funcs[statistic],axis=1)
		
		#store dataset value of threshold
		thresh = (series.describe(
							percentiles=np.arange(0,1,.01)
						)
						.loc['0%':]
						.iloc[threshold]
				)
		
		#determine direction of threshold
		#if greater than 50, greater than or equal to,
		#otherwise, less than
		if higher:
			return markers[series.where(series > thresh).notnull()]
		else:
			return markers[series.where(series <= thresh).notnull()]
		
	def markers_by_intrasample_variability(self,
		df,higher=False,threshold=80):
		"""
		Output markers below a certain variance threshold for 
		replicates within all samples

		Parameters
		----------
		df
			dataframe to subset
		higher {True,False}
			Indicates whether to subset by markers higher (True) than or less than/equal to (False) the threshold score
		threshold [0,100]
			Quantile for statistic threshold

		Output
		------
		subsetted markers

		"""

		#transpose
		df = df.T

		samp_rep_dict = self.sample_replicate_dictionary

		markers = pd.Index([])

		for samp in samp_rep_dict.keys():
			reps = samp_rep_dict.get(samp)
			sub = df.loc[np.asarray(reps)]
			markers = markers.append( \
				self.markers_by_statistic_threshold(df=sub.T,
					higher=False,
					statistic='variance',
					threshold=threshold))

		return markers.unique()

	#HELPER FUNCTIONS
	def make_df_samples(self,df,agg='mean'):
		"""
		From df_replicates, make df_samples by applying to the replicates for a sample the estimator indicated in the helper function

		Parameters
		----------
		df
			replicates dataframe to aggregate
		agg: {'mean','median', 'variance'}
			statistic to apply to replicate values for sample membership

		Output
		------
		subsetted dataset
		"""
		df_replicates = df.copy()
		dictionary = self.sample_replicate_dictionary

		# initializing empty, desired dataframe of samples by proteins
		df_samples = pd.DataFrame(index=self.markers,columns=self.samples)

		#for each sample, apply the helper function to df_replicates sub dataframe to make df_samples
		for sample in dictionary.keys():

			#dataframe of replicates for a given sample
			holder = df_replicates.loc[:,dictionary[sample]]

			#determine sample protein value from replicate values
			df = holder.agg(func=agg,axis=1)

			#put in df to larger df_samples
			df_samples.loc[:,sample] = df

		df_samples.columns.name = "sample"
		df_samples.index.name = self.marker_type

		return df_samples

class IntegratedCohort(object):
	"""
	Integrated patient cohorts object for shared patient clinical and proteomics data

	Dataframes, variables, and functions facilitating the processing, analysis, and integration of the integrated cohort data.
	
	Parameters
	----------
	cohort_objs
		dictionary of individual cohort objects to integrate through shared clinical and proteomic data

	Examples
	--------

	"""

	def __init__(self,cohort_objs={}):

		self.cohort_objs = cohort_objs
		self.df = None
		self.df_groups = None

	def sayHello(self):

		return "Hello!"

	def integrate_cohorts(self,treat_name="PGD",ref_name="NL",dataset_type='samples_hq',groups_type='sample_groups'):
		'''
		For created dataset integrating cohorts data
		
		Uses treat and ref attributes of each object to integrate the sample groups. Uses treat_name and ref_name to rename the treat and ref groups (in case in cohorts the treat and ref groups are named differently but should be integrated)
		
		Can set the object attribute type to use for joining e.g. the samples or replicates dataset
		
		Returns:
		--------
		
		integrated_samples_dataset:
			common proteins x all cohort samples
		integrated_sample_groups
			groups x all cohort samples
		response
			response_group x samples
		'''
		dfs = {}
		df_groups = {}
		for cohort,obj in self.cohort_objs.items():
			dfs[cohort] = getattr(obj,dataset_type)
			df = getattr(obj,groups_type)
			df_groups[cohort] = df

		self.df = multi_df_join(dfs).dropna()
		self.df_groups = multi_df_join(df_groups).dropna()

def make_cohorts_dict(names,file_dirs,replicates_files,samples_file,sample_groups_files,data_dirs,marker_types,references,treatments):
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
								data_dir=data_dirs[cohort],file_dir=file_dirs[cohort],
								replicates_file=replicates_files[cohort],
								samples_file=replicates_files[cohort],
								sample_groups_file=sample_groups_files[cohort],
								marker_type=marker_types[cohort],
								reference=references[cohort],treatment=treatments[cohort])
			objs[cohort] = obj
		return objs

if __name__=="__main__":
	"""
	Execute
	"""
	name =  'cohort_name'
	c = cohorts.Cohort(name)
