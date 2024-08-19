## Code for identifying top and bottom expressors of SLC7A11 in GTEx skeletal muscle data


# import modules
import numpy as np
import pandas as pd
import collections

import conorm 	# TMM-normalization

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_validate, StratifiedKFold
from sklearn.metrics import r2_score
from skrebate import TuRF, ReliefF
import scipy.stats as stats
from statsmodels.stats.multitest import fdrcorrection

# Define a funtion to format the age-brackets provided for GTEx samples to ge the minimum age
def format_age(sample_info):
	sample_info['Age'] = sample_info.loc[:,'Age Bracket'].str.split('-', n=1).str[0] 	# splits age range into min and max values, then return the minimum age value as a new column
	sample_info['Age'] = (sample_info['Age']).astype(int)
	return sample_info

# Given a subset of histology sample information, return a corresponding subset of gene expression data
def subset_expr(subset_ids, expr_df):

	uniq_ids = subset_ids.loc[:, 'Subject ID'].unique() 	# histology file contains sample information for all tissues, thus remove redundant data for our use

	subset_columns = [col for col in expr_df.columns if any(id_str in col for id_str in uniq_ids)]

	return expr_df[subset_columns]

# define a method to detect and remove outliers by Z-score before selecting ourtop and bottom expressors
def remove_outliers(expr_data):

	# Get Z-scores
	Z_expr = (expr_data - expr_data.mean()) / expr_data.std()

	# Return expression data with outliers removed, identified as subjects with Z-score > 3 or < -3
	return Z_expr[Z_expr.between(-3,3)]

# define a method to sort and subset expression data into top and bottom expressors of a spceified gene
def get_expressors(expr_df, sort_by, n_samples=25, clean_outliers=True):

	# get expression data for the gene we want to sort and select by, and removee outlier samples
	gene_expr = expr_df.loc[sort_by]
	if clean_outliers:
		gene_expr = remove_outliers(gene_expr)

	top_expr = gene_expr.nlargest(n_samples).index
	bottom_expr = gene_expr.nsmallest(n_samples).index

	return top_expr, bottom_expr

# Want to examine the gene expression distributions for each of our groups to observe any differences
# Given a wide format dataframe of genes x samples, convert to a longform DF and plot histogram of expression values
def plot_exprDistribution(expr_df, group_ids=False, group_names=False, group_colours=False):

	# if groups are specified, first split df into groups for plotting
	group_dfs=[]
	if group_ids:
		for ids in group_ids:
			temp_df = expr_df.loc[:, ids]
			temp_flat = temp_df.values.flatten()
			group_dfs.append(temp_flat)

		plot_df = pd.DataFrame(data=group_dfs)

		if group_names:
			plot_df.index = group_names

	# if no groups are specified just flatten the whole expression dataframe
	else:
		plot_df = pd.DataFrame(data=expr_df.values.flatten())

	plot_df.replace(-np.inf, np.inf, inplace=True)

	plot_df = plot_df + 0.0001
	if group_colours:
		hist = sns.histplot(data=plot_df.T, bins=50, log_scale=(True,True), element="step", palette=group_colours)
	else:
		hist = sns.histplot(data=plot_df.T, bins=50, log_scale=(True,True), element="step")

	return hist

# Define a function to create a split barplot for our expression groups through multiple axes to create broken y-axis
def split_bar(plot_data, plot_x, plot_y, plot_hue=False):
	fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex=True)

	bar1 = sns.barplot(data=plot_data, ax=ax1, x=plot_x, y=plot_y, hue=plot_hue, palette=['blue', 'orange'], dodge=False)
	bar2 = sns.barplot(data=plot_data, ax=ax2, x=plot_x, y=plot_y, hue=plot_hue, palette=['blue', 'orange'], dodge=False)

	# Adjust x-axis to only show certain samples names
	bar2.get_xaxis().set_visible(False)
	# bar2.xaxis.set_major_locator(ticker.MultipleLocator(3)) 	# show every 5th label
	plt.xticks(rotation=90)

	bar1.set_ylim(0.1,0.6) 	# adjust y-axis for each subplot
	bar2.set_ylim(0.00,0.04)

	# remove legends for each plot and create single axis outside
	bar1.get_legend().remove()
	bar2.get_legend().remove()

	bar2.legend(loc=(1.05,1), title='GTEx samples')

	# Adjust y-labels
	bar1.set_ylabel('')
	bar2.set_ylabel('')

	bar1.get_xaxis().set_visible(False) 	# x-axis is shared with lower plot

	return fig

# define a method for identifying significant genes and differential expression as volcano plot
def volcano(data, x, y, x_thresh=1, y_thresh=0.05, sig_palette=None):
	"""
	data = input gene expression dataframe (genes as rows)
	x = Column name containg fold change values
	y = Column name containing p-values
	x_thresh = significance threshold for fold change values
	y_thresh = significance threshold for p-values
	sig_palette = list of colour values for [negative DEGs, postive DEGs]
	"""

	# Want to identify the significant genes according to our x_thresh, y_thresh vals
	pos_genes =	data.where((data[y] < y_thresh) & (data[x] >= abs(x_thresh)))
	pos_genes.dropna(axis=0, how='all', inplace=True)
	neg_genes =	data.where((data[y] < y_thresh) & (data[x] <= -1*abs(x_thresh)))
	neg_genes.dropna(axis=0, how='all', inplace=True)
	print(f'pos_genes:\n{pos_genes}\nneg_genes:\n{neg_genes}')
	# define gene groups columns for colouring volcano plot points
	data['Differential expression'] = 'Non-significant'
	data.loc[pos_genes.index, 'Differential expression'] = 'Positive DEG'
	data.loc[neg_genes.index, 'Differential expression'] = 'Negative DEG'
	print(f'pos genes:\n{data.loc[pos_genes.index, "Differential expression"]}\n\nneg_genes:\n{data.loc[neg_genes.index, "Differential expression"]}')
	data['-log10 p-value'] = np.log10(data[y])*-1 	# format values for y-axis

	# Define colour palette for plot
	if sig_palette:
		plot_palette = {
			'Negative DEG': sig_palette[0], 
			'Non-significant': 'darkgrey',
			'Positive DEG': sig_palette[1]
			}
	else:
		plot_palette = {
			'Negative DEG': 'tomato', 
			'Non-significant': 'darkgrey',
			'Positive DEG': 'royalblue'
			}

	# # palette is coloured according to their appearance in the dataset
	# # Can guarentee that our ordering is correct by sorting the data values before plotting
	# sort_data = data.sort_values(by=x, ascending=True)

	# plot the data
	volc = sns.scatterplot(data=data, 
			x=x, 
			y='-log10 p-value',
			hue='Differential expression',
			alpha=0.70,
			linewidth=0,
			s=10,
			palette=plot_palette
			)

	# Set x_lim such that no extreme, erroneuous log2FC values are shown
	# Some |FC values| > ~8 appear as "extra tails" due to some very low (but non-zero) expression
	volc.set_xlim([-6,6])
	volc.set_ylim([-0.2,data['-log10 p-value'].max()+0.2])
	# Add lines corresponding to significance thresholds
	ymin, ymax = volc.get_ylim()
	xmin, xmax = volc.get_xlim()

	volc.vlines(x=[-1*abs(x_thresh),abs(x_thresh)], ymin=ymin, ymax=ymax, colors='grey', ls='--', lw=0.75)
	volc.hlines(y=-1*np.log10(y_thresh), xmin=xmin, xmax=xmax, colors='grey', ls='--', lw=0.75)

	plt.legend(bbox_to_anchor=(1.04,1), loc='upper left')

	return volc

# define method for computing VIP scores for PLS-DA
# from https://github.com/scikit-learn/scikit-learn/issues/7050
def VIP(model):
    t = model.x_scores_
    w = model.x_weights_ # replace with x_rotations_ if needed
    q = model.y_loadings_ 
    features_, _ = w.shape
    vip = np.zeros(shape=(features_,))
    inner_sum = np.diag(t.T @ t @ q.T @ q)
    SS_total = np.sum(inner_sum)
    vip = np.sqrt(features_*(w**2 @ inner_sum)/ SS_total)
    return vip

# define method to format inputs and run PLS-DA on data to extract important features
def run_plsda(data, classes, standardize=True, get_weights=True, get_vip=True):

	"""
	data = input gene expression data (n_genes, n_GTExSamples)
	classes = GTEx samples for categories we wish to discrimenate [class_0, class_1]
	standardize = whether or not the data need to be standardized before running PLS-DA (like PCA, PLS-DA inputs should be standardized)
	get_weights = Whether or not we want to return the feature weights importance instead of the PLSDA object 
	"""

	# Select our data such that we have columns = class 0, class 1 and generate a corresponding class vector
	x = data.loc[:, classes[0] + classes[1]]

	if standardize:
		scaler = StandardScaler()
		x = scaler.fit_transform(x.T) 	# n_samples, n_features

	y = np.array([[1]*len(classes[0]) + [0]*len(classes[1]), 	# create dummy 2,n array for membership in binary classes across samples
	[0]*len(classes[0]) + [1]*len(classes[1])]
	)
	print(f'input array x:\n{x}\n x_shape = {x.shape}')
	print(f'target labels array:\n{y}\ny_shape = {y.shape}')

	# Generate a model PLS-DA
	pls = PLSRegression(scale=False, n_components=2)

	# Cross-validate the model and return the mean Q^2 score over CVs to evaluate our model
	rng = np.random.default_rng(seed=42)
	idxs = rng.permutation(y[0].shape[0]) 	# get and shuffle indexs

	cv = cross_validate(pls, x[idxs], y[:,idxs].transpose(), scoring='r2', return_train_score=True, cv=5)
	print(f'Keys: {cv.keys()}')
	scores = [cv['train_score'], cv['test_score']]

	# Fit model to all of our data
	pls.fit(X=x, Y=y.transpose())

	if get_weights and get_vip:
		weights = pls.x_weights_
		vips = VIP(pls)
		print(f'vip shape: {vips.shape}')
		vips = vips.reshape(-1,1)
		print(f'new shape: {vips.shape}')
		return weights, scores, vips
	elif get_weights:
		return pls.x_weights_, scores

	elif get_vip:
		vips = VIP(pls)
		print(f'vip shape:\n{vips}')
		return scores, vips

	else:
		return pls.transform(X=x, Y=y.tranpose()), scores




# define a method to run skrebate's Relief methods for feature importance/selection
def run_relief(data, targets, discrete=False, multiround=False, multi_headers=None, evaluate=False):
	'''
	data (np.array): 	array of shape n_samples, n_features
	targets (np.array):	1D array of the targets for relief algorithm, can be discrete ints for classes
						or continuous for regression
	discrete (bool):	Sepcify type for 'targets', whether discrete (classes) or 
	multiround (bool):	Whether to implement single-round selection (ReliefF) or multi-round (TuRF). 
						For feature spaces > 10,000 multi-round is recommended
						(https://epistasislab.github.io/scikit-rebate/using/#general-usage-guidelines)
	'''

	# # Select our data such that we have columns = class 0, class 1 and generate a corresponding class vector
	# x = data.loc[:, classes[0] + classes[1]]

	# if standardize:
	# 	scaler = StandardScaler()
	# 	x = scaler.fit_transform(x.T) 	# n_samples, n_features

	# y = np.array([[1]*len(classes[0]) + [0]*len(classes[1]), 	# create dummy 2,n array for membership in binary classes across samples
	# [0]*len(classes[0]) + [1]*len(classes[1])]
	# )

	print(f'Standardizing')
	print(f'x shape: {data.shape}, y sahpe: {targets.shape}')
	scaler = StandardScaler()
	data = scaler.fit_transform(data)  	# n_samples, n_features
	# targets = scaler.fit_transform(targets)
	targets = (targets - targets.mean())/targets.std()
	print(f'After scaling:\nx shape: {data.shape}, y sahpe: {targets.shape}')
	if multiround:
		fs = TuRF(core_algorithm="ReliefF", n_features_to_select=10, pct=0.5, verbose=True, n_neighbors=20)  	# NOTES: 1) pct > 0.5 appears to cause error (see https://github.com/EpistasisLab/scikit-rebate/issues/54)
																								# 		 2) If number of features are odd, pct splits will be uneven and result in missing idxs and error
		fs.fit(data, targets, multi_headers)

	else:
		# fs = ReliefF()
		# fs.fit(data, targets)

		pass

	scores = fs.feature_importances_

	if evaluate:
		print(f'targets shape: {targets.shape}')
		# Cross-validate the model and return the mean Q^2 score over CVs to evaluate our model
		rng = np.random.default_rng(seed=42)
		idxs = rng.permutation(targets[0].shape) 	# get and shuffle indexs

		fs = ReliefF()
		cv = cross_validate(fs, data[idxs], targets[idxs].transpose(), scoring='r2', return_train_score=True, cv=5)
		print(f'Keys: {cv.keys()}')
		evals = [cv['train_score'], cv['test_score']]

		return scores, evals


		### 	!!!!!!!!!!!!!!!!!!!!!!!!!!!
		###
		###		NOTE: AttributeError: 'ReliefF' object has no attribute 'predict', NO CV POSSIBLE
		###
		###		!!!!!!!!!!!!!!!!!!!!!!!!!!!

	else:
		return scores



# define recursive function for flattening a list which contains lists and int values such that the return is a list of all elements
def flatten(xs):
	"""
	xs = a list of lists and values
	"""
	if isinstance(xs, collections.Iterable):
		return [a for x in xs for a in flatten(x)]
	else:
		return [xs]


# recursive function to call in get_nearest to get original values from clusters in a linkage matrix
def get_child(Z, i, n, idxs=False):

	"""
	Z = a scipy (or sklearn?) linkage matrix where Z[i] = [node_j, node_k, distance, n_values]
	i = index where we need to get child leaves
	n = number of original values, will inform what is a new cluster vs. original value
	idxs = Optional, list of linkage matrix idxs to append child idxs to
	"""
	# Get our index values to return
	jeff = 0
	val_0 = Z[int(i)][0]
	val_1 = Z[int(i)][1]
	idx_0 = int(i)
	idx_1 = int(i)

	# Check if these values are original values or clusters (where val > n)
	if val_0 > n:
		# print(f'!!get_child()!! val_0 > n')
		idx_0, val_0 = get_child(Z, val_0-n, n, idxs)
	if val_1 > n:
		# print(f'!!get_child()!! val_1 > n')
		idx_1, val_1 = get_child(Z, val_1-n, n, idxs)
	else:
		# print(f'!!get_child()!! all vals < n')
		jeff = 1

	# return the indexes obtained over recursively getting clusters until original values
	if idxs:
		idxs = flatten([idx_0, idx_1])
		# print(f'get_child() return, i={int(i)}:\n\tl={Z[int(i)]}\n\tReturns:\n\t\tidxs: {idxs}\n\t\tvals: {[val_0,val_1]} ')
		return idxs, [val_0, val_1]

	else:
		return [val_0, val_1]


# Define function to get the nearest neighbours of a given observation in a hierarchical clustering linkage map
def get_nearest(Z, node_index, max_steps=False, max_dist=False, get_links=False):

	"""
	Z = a scipy (or sklearn?) linkage matrix where Z[i] = [node_j, node_k, distance, n_values]
		is the new linkage formed at step_i in the clustering algorithm, where a node may be an original value or a new cluster if node_index > n (original samples) and node index = n+i
	node_index = value index of either the original value we are interested in, or a new cluster formed by the algorithm if node_index > n
	max_steps = thresholding parameter for the number of clustering steps to extend our nearest neighbour search for
	max_dist = thresholding parameter for a distance cut-off to threshold the number of clusters we obtain
	sub_links = whether to get and return the sub linkage matrix, 
	"""

	n = np.shape(Z)[0] + 1 	# number of original observations/values, corresponds to rows + 1 in Z
	links_idx = [] 	# list of all linkage matrix indexes for generating a sub matrix if get_links=True
	steps = 0 	# Thresholding checks
	last_i = 0

	# Find the first clustering of our node of interest and get all original values in this clustering (using get_child())
	for i, l in enumerate(Z):
		# First, iterate over the matrix and work "up" the linkage matrix while checking for our thresholding parameters
		if l[0] == node_index or l[1] == node_index:
			links_idx.append(i)

			if max_steps and steps < max_steps: 	# Check that we have not reached max steps
				steps += 1
				node_index = i + n

			elif max_dist and l[2] < max_dist: 	# Check that we have not surpassed max distance
				last_i = i
				node_index = i + n

			else:
				if max_dist:
					l = Z[last_i] 	# For max_dist threshold, we check that our cluster distance is less than thershold, when False we must proceed with the last cluster under the distance threshold
					del links_idx[-1]
				# once steps = max_steps we can work back through the linkage matrix to get all the original values contained within clusters
				if l[0] > n:
					# print(f'[l[0]] > n')
					idxs_0, val_0 = get_child(Z, l[0]-n, n, links_idx)
				else:
					val_0 = l[0]

				if l[1] > n:
					# print(f'[l[1]] > n')
					idxs_1, val_1 = get_child(Z, l[1]-n, n, links_idx)
				else:
					val_1 = l[1]

	# convert val_0 and val_1 from lists of lists to flattened list for easy indexing
	vals = [x for xs in [val_0, val_1] for x in xs]

	# now can return the final list of values
	if get_links:
		idxs = [x for xs in [links_idx, idxs_0, idxs_1] for x in xs]
		idxs = sorted(list(set(idxs)))
		return Z[idxs], vals
	else:
		return vals

# define function for applying PCA to our data
def fit_pca(expr_data):

	# define PCA params and fit data
	pca_alg = PCA(n_components=2)
	princip_componenets = pca_alg.fit_transform(expr_data) 	# fit_transform takes rows as samples and cols as features, so this will result in gene-wise PCA

	var = pca_alg.explained_variance_ratio_

	col_names = [f'Principal component 1 ({var[0]*100:.1f}%)', f'Principal component 2 ({var[1]*100:.1f}%)']

	PC_df = pd.DataFrame(data=princip_componenets, index=expr_data.index, columns=col_names)

	return PC_df

# define a wilcoxon method which can be applied between our groups across genes to determine significant differences
# Test of independence is wilcoxon ranked sums test
def gene_rank_sum(geneExpr, group_a, group_b):
	geneExpr_a = geneExpr.loc[group_a]
	geneExpr_b = geneExpr.loc[group_b]

	# Apply wilcoxon test
	stat, pval = stats.ranksums(x=geneExpr_a, y=geneExpr_b, alternative='two-sided')

	return pval

# define a method for determining FC between our groups (Fold change as group_a/group_b)
# Ideally should be able to determine multiple FCs (log10, log2, median, mean, etc.)
def fold_change(geneExpr, group_a, group_b, method='median', log2=True):

	geneExpr_a = geneExpr.loc[:, group_a] 	# get gene expression for each group
	geneExpr_b = geneExpr.loc[:, group_b]

	# Get representative gene expression values for each group based on {method}
	if method == 'median':
		comp_a = geneExpr_a.median(axis=1)
		comp_b = geneExpr_b.median(axis=1)

	elif method == 'mean':
		comp_a = geneExpr_a.mean(axis=1)
		comp_b = geneExpr_b.mean(axis=1)

	# Caclulate fold change for each group
	gene_fc = (comp_a + 0.000001) / (comp_b + 0.000001) 	# Add small constant to avoid 0/0

	# convert to log FC
	if log2 == True:
		logFC_expr = pd.Series(data=np.log2(gene_fc), index=geneExpr.index, name='Log2 fold change')
	else:
		logFC_expr = pd.Series(data=np.log(gene_fc), index=geneExpr.index, name='Log2 fold change')

	return logFC_expr


# # get top and bottom n expressors of SLC7A11 for each of our groups
# slc7a11_id = 'ENSG00000151012.13'
# old_exprDF = old_exprDF.loc[:, old_exprDF.loc[slc7a11_id] !=0] 	# remove any 0 expressors before getting top and bottom expressors
# young_exprDF = young_exprDF.loc[:, young_exprDF.loc[slc7a11_id] !=0]

# old_top, old_bottom = get_expressors(expr_df=old_exprDF, n_samples=25, sort_by=slc7a11_id)
# young_top, young_bottom = get_expressors(expr_df=young_exprDF, n_samples=25, sort_by=slc7a11_id)

# ## Create histograms to understand the distribution of gene expression in our different groups
# expr_plt = plot_exprDistribution(expr_df=tmm_exprDF, group_ids=[old_top, old_bottom], group_names=['Top 25 Old', 'Bottom 25 Old'], group_colours=['blue', 'orange'])
# expr_plt.set_xlabel('mRNA Expression (TPM, TMM normalized)')
# fig = expr_plt.get_figure()
# fig.savefig('expr_hist_loglog_TMM_old.png')
# plt.close(fig)

# expr_plt = plot_exprDistribution(expr_df=tmm_exprDF, group_ids=[young_top, young_bottom], group_names=['Top 25 Young', 'Bottom 25 Young'], group_colours=['purple', 'gold'])
# expr_plt.set_xlabel('mRNA Expression (TPM, TMM normalized)')
# fig = expr_plt.get_figure()
# fig.savefig('expr_hist_loglog_TMM_young.png')
# plt.close(fig)

# expr_plt = plot_exprDistribution(expr_df=tmm_exprDF, group_ids=[old_top, old_bottom, young_top, young_bottom], group_names=['Top 25 Old', 'Bottom 25 Old', 'Top 25 Young', 'Bottom 25 Young'], group_colours=['blue', 'orange', 'purple', 'gold'])
# expr_plt.set_xlabel('mRNA Expression (TPM, TMM normalized)')
# fig = expr_plt.get_figure()
# fig.savefig('expr_hist_loglog_TMM_all.png')
# plt.close(fig)


# plot_data = tmm_exprDF.loc[slc7a11_id, old_top.tolist() + old_bottom.tolist()[::-1]] 	# get expression values for only the top and bottom IDs (NOTE: reverse bottom list to be ordered largest-> smallest)

# # plot these expression levels as a bar plot
# plot_df = plot_data.to_frame()

# plot_df['GTEx samples'] = 'Top 25 Old'
# plot_df.loc[old_bottom.tolist(),'GTEx samples'] = 'Bottom 25 Old'
# plot_df['GTEx sample ID'] = plot_df.index
# plot_df.to_csv(rf'.\data\GTEx\oldHighLow_xctExpressors_GTExv8.csv')
# # exit()

# slc7a11_expr = split_bar(plot_data=plot_df, plot_x='GTEx sample ID', plot_y='ENSG00000151012.13', plot_hue='GTEx samples')
# slc7a11_expr.text(0.0, 0.55, 'SLC7A11 muscle mRNA (TPM, TMM normalized)', va='center', rotation='vertical') 	# Add y-axis text across both subplots
# slc7a11_expr.savefig(f'tmmSLC7A11_highLowGTEx_old_Broken.png', bbox_inches='tight')

# plot_data = tmm_exprDF.loc[slc7a11_id, young_top.tolist() + young_bottom.tolist()[::-1]] 	# get expression values for only the top and bottom IDs (NOTE: reverse bottom list to be ordered largest-> smallest)

# # plot these expression levels as a bar plot
# plot_df = plot_data.to_frame()

# plot_df['GTEx samples'] = 'Top 25 Young'
# plot_df.loc[young_bottom.tolist(),'GTEx samples'] = 'Bottom 25 Young'
# plot_df['GTEx sample ID'] = plot_df.index
# plot_df.to_csv(rf'.\data\GTEx\youngHighLow_xctExpressors_GTExv8.csv')
# # # exit()

# slc7a11_expr = split_bar(plot_data=plot_df, plot_x='GTEx sample ID', plot_y='ENSG00000151012.13', plot_hue='GTEx samples')
# slc7a11_expr.text(0.0, 0.55, 'SLC7A11 muscle mRNA (TPM, TMM normalized)', va='center', rotation='vertical') 	# Add y-axis text across both subplots
# slc7a11_expr.savefig(f'tmmSLC7A11_highLowGTEx_young_Broken.png', bbox_inches='tight')

