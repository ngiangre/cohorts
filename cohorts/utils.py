import numpy as np
import pandas as pd
import logging
from functools import reduce
import scipy.stats as sc

def allq(X):
	"""
	Helper function used in manual_feature_extraction

	Parameters
	----------
	X
		dataframe to subset

	Outputs
	-------
	Subsetted dataframe with all proteins quantified
	"""
	applied = X.apply(lambda x: np.sum(x>0),axis=1)
	bools = np.where(applied == len(X.columns))
	return X.iloc[bools].index.tolist()

def allnotq(X):
	"""
	Helper function used in manual_feature_extraction

	Parameters
	----------
	X
		dataframe to subset

	Outputs
	-------
	Subsetted dataframe with no proteins quantified
	"""
	applied = X.apply(lambda x: np.sum(x),axis=1)
	bools = np.where(applied == 0)
	return X.iloc[bools].index.tolist()

def mixed(X):
	"""
	Helper function used in manual_feature_extraction

	Parameters
	----------
	X
		dataframe to subset

	Outputs
	-------
	Subsetted dataframe with either 1 up to all proteins quantified
	"""
	applied = X.apply(lambda x: np.sum(x>0),axis=1)
	condition1 = (applied < len(X.columns))
	condition2 = (applied > 0)
	bools = np.where( condition1 & condition2 )
	return X.iloc[bools].index.tolist()

def get_uniprot_table(file="/Users/npg2108/Research/Projects/exosome_pgf/data/uniprot-all_20171124.tab.gz"):
	"""
	Uploading uniprot database flat file
	"""
	tab = pd.read_csv(file,delimiter="	",index_col=0)

	return tab

def get_real(series, agg='mean'):
	"""
	Helper function to get array mean or return NaN for series of NaNs/floats. 
	Used to get sample value of protein for replicates in df_replicates
	"""

	func = {'mean' : np.mean,'median': np.median,'variance' : np.var}

	# get logical list of whether the protein value is NaN
	inds = [np.isnan(x) for x in series]
	
	#how many are NaN?
	summ = np.sum(inds)
	
	#give NaN/float for protein value depending on number of NaNs in series
	#returns NaN if all NaN values, and function value if atleasst one is non-NaN
	if summ==len(inds):
		return 0
	else:
		#negativing-list of booleans for indices where non-NaN
		reals = [not i for i in inds]
		stat = func[agg](series[ reals ])
	return stat

def rank_normalize(X):
	"""
	Helper function to normalize values in sample by normalized rank, so that protein values within a sample are [0,1] and indicate relative expression in that sample. Hopefully this will allow better comparisons between samples in a cohort.
	"""

	df_ranked = X.apply(lambda x : sc.rankdata(x,method='dense'))

	df_ranked_normalized = df_ranked.apply(lambda x : (x / np.max(x)))

	return df_ranked_normalized 

def make_comparable(dfs,rank=True):
	"""
	Make list of cohort datasets (proteins x reps/samps) comparable i.e. make sure the 
	proteins are the same between the datasets and normalize.
	"""

	#get list of sets of proteins in each cohort dataset
	if type(dfs) == dict:
		prot_set_lst = [set(x.index) for i,x in dfs.items()]
		dfs = [x for i,x in dfs.items()]
	else:
		prot_set_lst = [set(x.index) for x in dfs]

	#get protein intersection
	com_prots = multiple_set_intersection(*prot_set_lst)

	#iterate through dfs and normalize and append to new list
	dfs_new = []
	for df in dfs:
		if rank:
			dfs_new.append(rank_normalize(df.loc[com_prots]))
		else:
			dfs_new.append(df.loc[com_prots])

	return dfs_new 

def multiple_set_intersection(*sets):
    """
    Return multiple set intersection.
    from https://stackoverflow.com/questions/2541752/best-way-to-find-the-intersection-of-multiple-sets
    """
    try:
        return set.intersection(*sets)
    except TypeError: # this is Python < 2.6 or no arguments
        pass

    try: a_set= sets[0]
    except IndexError: # no arguments
        return set() # return empty set

    return reduce(a_set.intersection, sets[1:])

def quantileNormalize(df_input):
	"""
	From https://github.com/ShawnLYU/Quantile_Normalize/blob/master/quantile_norm.py
	Quantile normalization method of a dataset
	"""
	df = df_input.copy()
	#compute rank
	dic = {}
	for col in df:
		dic.update({col : sorted(df[col])})
	sorted_df = pd.DataFrame(dic)
	rank = sorted_df.mean(axis = 1).tolist()
	#sort
	for col in df:
		t = np.searchsorted(np.sort(df[col]), df[col])
		df[col] = [rank[i] for i in t]
	return df

	#from here
#https://gist.github.com/mpschr/5db20df78c034654f030

def multi_df_join(df_dict):
	"""
	Takes a list of dataframes with the same primary key (id) and merges the columns, using
	the dictionary keys as suffixes in case there are conflicts with the column names.
	Parameters
	----------
	df_dict: A dictionary where keys are strings (collision suffixes) and values are DataFrames
	
	Output
	------
	The merged data frames
	"""

	# check if there are colliding column names
	all_columns = []
	for df in df_dict.values():
		all_columns = all_columns + list(df.columns)
		colliding_columns = set([x for x in all_columns if all_columns.count(x) > 1])

	df = None
	for suffix, input_df in df_dict.items():

		# rename colliding columns
		renamer_dict = {}
		for col in list(input_df.columns):
			if col in colliding_columns:
				renamer_dict[col] = col + suffix
		input_df.rename(columns=renamer_dict, inplace=True)

		# join columns
		if df is None:
			df = input_df
		else:
			df = df.join(input_df, how='outer')
		logging.debug("Shape is currently: {}".format(df.shape))

	return df

def treat_ref_color_map(df,labels,groups,palette='hls'):
	'''
	Create color dictionary for reference and treatment samples for use in plotting

	Parameters
	----------
	df
		group labels from sample or replicate groups
	labels
		array of ref and treat group label
	groups
		label from dataframe in which we want to attribute colors to
	Order by the labels

	Output
	------
	Samples to rgb pandas Series
	'''

	#create array of group names for correct sample membership
	cond = [] * df.shape[0]
	for g in groups:
		for n in df.loc[g].ravel():
			if n == 1:
				cond.append(g)

	#create a color palette with the same number of colors as unique groups-right now only 2
	network_pal = sns.color_palette(palette,n_colors=len(groups))

	#create a dictionary where the key is the group and the values are the colors from the palette we just created
	network_lut = dict(zip(groups, network_pal))

	#convert array of group names to sample membership to a pandas series for mapping
	networks = pd.Series(cond,index=labels)

	#map the colors to the series. Now we have a list of colors the same length as our dataframe, where unique values are mapped to the same color
	network_colors = networks.map(network_lut)

	#change index class to string
	network_colors.index = network_colors.index.astype('str')

	return network_colors


