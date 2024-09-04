"""
Luke Kennedy - Fong-McMaster et al. (2024)

Code for identifying differentially expressed genes between pre-specified groups of OPA1 expressors as reported in Fong-McMaster et al. (2024)
	Given a csv file containing top and bottom expressors as determined by `get_expression_groups.py`, perform statistical differential 
	expression analysis, PLS-DA, and Relief analysis. Output results for each and accompanying figures
"""


# import modules

import get_GTEx  	# helper functions

import numpy as np
import pandas as pd
import expression_analysis 		# Helper functions
from statsmodels.stats.multitest import fdrcorrection

# Plotting
import seaborn as sns
import matplotlib.pyplot as plt

# define main function for our analyses
#	Subsets whole GTEx data to just get left ventrical data -> output leftVent_{GTEX}.csv 	(mRNA TPM values)
#	From this, select sample IDs in top n or top n% expression of OPA1. Sex-controlled, ages 20-59
# 	Initial analysis includes a comparison of pathology/hardy scales between these two groups
def main():

	'''
	Arguments
	TOP_EXPR  (str): 				File path for our identified top OPA1 samples gene expression
	BOTTOM_EXPR  (str): 			File path for our identified bottom OPA1 samples gene expression

	GENE_ID (str):					Gene ID to subset our data by
 	COLOURS (list, strs):			Colours for plotting as hex codes
	''' 

	TOP_EXPR = r'\data\top_LV_OPA1_RNAseq.csv' 	# !!PLACEHOLDER!! Replace 'data' with your directory
	BOTTOM_EXPR = r'\data\bottom_LV_OPA1_RNAseq.csv' 	# !!PLACEHOLDER!! Replace 'data' with your directory
	GENE_ID = 'ENSG00000198836.8' 	# OPA1
	COLOURS = ['#F9A51B', '#3A54A4'] 	# Orange and blue


	# Get our expression data for each group, generated previously by `get_expression_groups.py`
	top_df = pd.read_csv(TOP_EXPR, index_col=0)
	bottom_df = pd.read_csv(BOTTOM_EXPR, index_col=0)

	'''
	To further compare our groups, get differentially expressed genes between our groups
	Define differential expression as median log2FC > |1| and BH-adjusted Wilcoxon p-values < 0.05
	'''
	# Get differentially expressed genes between high/low groups
	groups_exprDF = pd.concat([top_df, bottom_df], axis=1)

	p_vals = groups_exprDF.apply(lambda x: expression_analysis.gene_rank_sum(geneExpr=x, group_a=top_df.columns.tolist(), group_b=bottom_df.columns.tolist()), axis=1)

	groups_exprDF['wilcoxon p-values'] = p_vals
	groups_exprDF['BH adjusted p-values'] = fdrcorrection(p_vals.values, method='indep', alpha=0.05)[1] 	# get adjust p-values

	fc = expression_analysis.fold_change(geneExpr=groups_exprDF, group_a=top_df.columns.tolist(), group_b=bottom_df.columns.tolist(), method='median', log2=True)
	groups_exprDF.loc[:,'Log2 fold change'] = fc
	groups_exprDF.to_csv(rf'LV_GTEx_OPA_DEG.csv')

	# Now get significantly differentially erxpressed genes with FC > |1.5|
	sigDEGs = groups_exprDF.where((groups_exprDF['BH adjusted p-values'] < 0.05) & ((groups_exprDF['Log2 fold change'] <= -1) | (groups_exprDF['Log2 fold change'] >= 1)))
	print(f'sigDEGs shape: {sigDEGs.shape}')
	sigDEGs.dropna(axis=0, how='all', inplace=True)
	print(f'sigDEGs non-zero shape: {sigDEGs.shape}')
	sigDEGs.to_csv(rf'LV_GTEx_OPA_significant_DEG.csv')


	## To visualize our differential expression, we can create a volcano plot
	## Though this is simply a scatterplot, we need to do some special formatting so our results are clear
	## Use the volcano plot helper function to do some formatting

	DE_plot = expression_analysis.volcano(data=groups_exprDF, 
			x='Log2 fold change', 
			y='BH adjusted p-values', 
			x_thresh=1.5, 	# specific our log2 FC threshold (absolute)
			y_thresh=0.01, 	# specify our FDR threshold
			sig_palette=COLOURS
			)

	fig = DE_plot.get_figure()
	fig.savefig(f'OPA1_quintiles_volcano.png', bbox_inches='tight', dpi=600)
	plt.close()

	'''
	Another method for determining relevant genes is via PLS-DA.
	Here, we define two classes of samples (high/low expressors) which the model aims to classify
	We extract the most important feautres (genes) for this seperation as most import features.
	'''

	weights, scores, vips = expression_analysis.run_plsda(data=groups_exprDF, 
		classes=[top_df.columns.tolist(), bottom_df.columns.tolist()],  	# NOTE: order of groups determines class 1 & class 0; i.e., flips the signs
		standardize=True, 
		get_weights=True,
		get_vip=True
		)

	mean_scores = [np.mean(metrics) for metrics in scores]
	print(f'mean score: {mean_scores}')
	joint = np.concatenate((weights,vips),axis=1)
	weights_df = pd.DataFrame(data=joint, index=groups_exprDF.index, columns=['Weights_1','Weights_2','VIP'])

	## We can apply an iterative approach to the important features selection, where features are dropped
	## by VIP score and PLS-DA is applied with this subset to observe changes in performance
	## When performance drops, that is the signal that minimal  features have been reached
	weights_df.drop(labels=GENE_ID,axis=0, inplace=True) 	# Remove OPA-1 to observe relevance of other genes
	max_vip = weights_df['VIP'].max()
	max_int = round(max_vip,1)
	VIP_thresholds = np.arange(0.9, max_int,0.1)[:-1] 	# Want at least 2 features or else cannot compute 2 omponents
	r2_scores, q2_scores = [], []
	for threshold in VIP_thresholds:
		subset_genes = weights_df[weights_df.loc[:,'VIP'] > threshold].index 	# get genes where VIP > threshold

		sub_exprDF = groups_exprDF.loc[subset_genes]
		# Run PLS-DA with subset
		_, cv_scores = expression_analysis.run_plsda(data=sub_exprDF,
				classes=[top_df.columns.tolist(), bottom_df.columns.tolist()],
				standardize=True,
				get_weights=True,
				get_vip=False
				)
		mean_scores = [np.mean(score) for score in cv_scores]
		r2_scores.append(mean_scores[0])
		q2_scores.append(mean_scores[1])


	# Create a plot to show performance over thresholds
	print(f'r2: {r2_scores}\nq2: {q2_scores}\nvips: {VIP_thresholds}')
	plot_df = pd.DataFrame(data=[VIP_thresholds,q2_scores], index=['VIP score threshold', 'mean Q2'])
	score_plot = sns.lineplot(data=plot_df.T,
		x='VIP score threshold',
		y='mean Q2',
		marker='o',
		color=COLOURS[1]
		)
	plt.ylim(min(q2_scores)-0.005, max(q2_scores)+0.005) 	# Add some buffer around our range

	# Add an indicator at our max point for identifying the best VIP cutoff
	# x index for max y values for stim and cue
	max_q2 = max(q2_scores)
	max_vip = VIP_thresholds[q2_scores.index(max_q2)] 	# get position of max q2 in list to get corresponding VIP
	print(f'max_vip: {max_vip}')

	if isinstance(max_vip, np.ndarray):
		max_vip = max_vip[1] 	# If multiple VIP scores with same performance, select the highest threshold

	# y min and max
	ymin, ymax = score_plot.get_ylim()

	# vertical lines
	score_plot.vlines(x=max_vip, ymin=ymin, ymax=max_q2, colors=COLOURS[1], ls='--', lw=1)


	plt.title('PLS-DA performance with different feature subsets')
	fig = score_plot.get_figure()
	fig.savefig('LV_OPA-PLSDA_VIP_line.png', dpi=400)

	weights_df.to_csv(rf'OPA1_DE_plsda_weights.csv')

	'''
	Similarily, we can use a Relief-based method for feature importance.
	Here, we can treat our data as a regression tast where we want to identify the most important
	genes for predicting OPA1 expression in our samples from top/bottom expressors
	'''

	target_expr = groups_exprDF.loc[GENE_ID].values 	# get OPA expression values for regression targets
	other_expr = groups_exprDF.drop(labels=GENE_ID, axis=0) 	# All other values as features
	gene_ids = other_expr.index.tolist()

	relief_scores = expression_analysis.run_relief(data=other_expr.values.transpose(), 
		targets=target_expr, 
		discrete=False,
		multiround=True,
		multi_headers=gene_ids
		)

	for feature_name, feature_score in zip(gene_ids, relief_scores):
		print(f'{feature_name}:    {feature_score}')

	relief_df = pd.Series(data=relief_scores, index=gene_ids, name='Relief Scores')
	relief_df.sort_values(ascending=False, inplace=True)
	relief_df.to_csv(rf'relief_scores_k20_scaled.csv')
main()