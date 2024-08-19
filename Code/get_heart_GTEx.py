"""
Luke Kennedy - Fong-McMaster et al. (2024)

Code for accessing and analyzing bulk GTEx RNAseq data for evaluation of OPA1 expression as reported in Fong-McMaster et al. (2024)
	Given a gene expression data file (bulk or tissue-specific), this code selects top and bottom expressing samples for a gene of interest.
	Certain selection criteria can be applied through the argument variables ("AGE_CUTOFF", etc.) and 
"""


# import modules

import get_GTEx

import numpy as np
import pandas as pd
import expression_analysis 		# Helper functions
import conorm 		# TMM-normalization
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import mannwhitneyu, ttest_ind

# Plotting
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import colorsys

# define main function for our analyses
#	Subsets whole GTEx data to just get left ventrical data -> output leftVent_{GTEX}.csv 	(mRNA TPM values)
#	From this, select sample IDs in top n or top n% expression of OPA1. Sex-controlled, ages 20-59
# 	Initial analysis includes a comparison of pathology/hardy scales between these two groups
def main():

	'''
	Arguments
	BULK_GTEx (str): 				File path for our bulk GTEx RNAseq file
	GTEx_INFO (str):				File path for our GTEx info file
	GTEx_TISSUE (str): 				Tissue ID to subset our data by

	GENE_ID (str):					Gene ID to subset our data by
	AGE_CUTOFF (int):				Maximum age we want in our data
	SEX_CONTROL (bool):				If we want our results to be balanced by sex
	CND_CONTROL (dict):				Another condition to control for, where dict key is the sample feature and dict value is the feature value
	N_CUTOFF (int, opt): 			Number of top/bottom expressor samples to return
 	COLOURS (list, strs):			Colours for plotting as hex codes
	''' 

	BULK_GTEx = r'D:\Claire_GTEx\data\GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'
	GTEx_INFO = r'GTEx_histologyNotes.csv'
	GTEx_TISSUE = 'Heart - Left Ventricle'

	GENE_ID = 'ENSG00000198836.8' 	# OPA1
	AGE_CUTOFF = 79
	CND_CONTROL = {'Hardy Scale': 'Ventilator case'}
	SEX_CONTROL = True
	N_CUTOFF = 24 		# 24 = opa 10% (x2) after outliers, vents only, all ages
	COLOURS = ['#F9A51B', '#3A54A4'] 	# Orange and blue

	# temp = get_GTEx.get_subjectExpr(gtex_expr=BULK_GTEx, gtex_samples=GTEx_INFO, gtex_tissues=GTEx_TISSUE)

	# # Save the file for further downstream analyses
	# print(f'temp df for just left ventricle data:\n{temp}\nsahpe:{temp.shape}')
	# temp.to_csv(rf'outputs/leftVent_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.csv')

	# # Save the file for further downstream analyses
	# print(f'temp df for just left ventricle data:\n{temp}\nsahpe:{temp.shape}')
	# temp.to_csv(rf'outputs/leftVent_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.csv')

	temp = pd.read_csv(rf'outputs/leftVent_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.csv', index_col=1)
	print(f'initial expression df:\n{temp}')

	'''
	Data may be skewed across samples, so reads for some samples are, in general, higher or lower.
	Normalize the TPM values usign TMM, then determine top and bottom IDs
	See if these IDs are different compared to non-standardized values
	'''

	temp = temp.iloc[:,2:] 	# remove metadata
	temp = temp.loc[(temp != 0).any(axis=1)] 	# remove any completely unexpressed genes

	expr_df = conorm.tmm(temp) 	# normalize expression data

	# Optionally, want to subset the range of sample ages
	if AGE_CUTOFF:
		sample_df = pd.read_csv(GTEx_INFO)
		format_sampleDF = expression_analysis.format_age(sample_info=sample_df) 	# Re-factor age brackets to minimum ages

		age_ids = format_sampleDF.loc[format_sampleDF.loc[:,'Age'] < AGE_CUTOFF]

		expr_df = expression_analysis.subset_expr(subset_ids=age_ids, expr_df=expr_df)

		print(f'Age-subsetted expr_df:\n{expr_df}')

	# Optionally, want to control for another feature in the data
	# E.g., cause of death which can be approximated by Hardy scale
	if CND_CONTROL:
		for key, values in CND_CONTROL.items():
			print(f'CND_CONTROL True:\nKey: {key}, values: {values}')
			cnd_ids = format_sampleDF.loc[format_sampleDF.loc[:, key] == values]

			print(f'cnd_ids:\n{cnd_ids}')
			expr_df = expression_analysis.subset_expr(subset_ids=cnd_ids, expr_df=expr_df)

		print(f'Condition-controlled dataset (subsetted on {CND_CONTROL}):\n{expr_df}')

	# Identify top and bottom expressing groups of our gene of interest

	# get top and bottom n expressors of GENE_ID
	gene_exprDF = expr_df.loc[:, expr_df.loc[GENE_ID] !=0] 	# remove any 0 expressors before getting top and bottom expressors

	# remove outliers
	temp = gene_exprDF.loc[GENE_ID]
	noOut_exprDF = expression_analysis.remove_outliers(temp) 	# Remove outliers (|Z-score| > 3)
	gene_exprDF = gene_exprDF.loc[:,noOut_exprDF.index] 		# subset the full exprdf to remove outliers
	print(f'gene_exprDF w/ no outliers:\n{gene_exprDF}')
	print(f'columns:\n{gene_exprDF.columns}')
	for ids in gene_exprDF.columns:
		print(ids)

	if SEX_CONTROL:
		male_ids = sample_df.loc[sample_df.loc[:,'Sex'] == 'male']
		female_ids = sample_df.loc[sample_df.loc[:,'Sex'] == 'female']

		male_exprDF = expression_analysis.subset_expr(subset_ids=male_ids, expr_df=gene_exprDF)
		female_exprDF = expression_analysis.subset_expr(subset_ids=female_ids, expr_df=gene_exprDF)

		print(f'Shapes of split data:\nMale={male_exprDF.shape}, female={female_exprDF.shape}, orig={gene_exprDF.shape}')

		top_male, bottom_male = expression_analysis.get_expressors(expr_df=male_exprDF, n_samples=N_CUTOFF, sort_by=GENE_ID, clean_outliers=False)
		top_female, bottom_female = expression_analysis.get_expressors(expr_df=female_exprDF, n_samples=N_CUTOFF, sort_by=GENE_ID, clean_outliers=False)

		top_maleDF = gene_exprDF.loc[:,top_male]
		bottom_maleDF = gene_exprDF.loc[:,bottom_male]
		top_femaleDF = gene_exprDF.loc[:,top_female]
		bottom_femaleDF = gene_exprDF.loc[:,bottom_female]
		print(f'Summary stats:\nTop male:\n{top_maleDF.describe()}\nBottom male:\n{bottom_maleDF.describe()}')
		print(f'Top female:\n{top_femaleDF.describe()}\nBottom female:\n{bottom_femaleDF.describe()}')

		# Combine groups
		top_df = pd.concat([top_maleDF,top_femaleDF],axis=1)
		bottom_df = pd.concat([bottom_maleDF,bottom_femaleDF],axis=1)

		print(f'summary satistics for combined groups:\nTop:\n{top_df.describe()}')
		print(f'Bottom:\n{bottom_df.describe()}')
		# exit()
		# Plot bar plots of our high/low expressors
		plot_data = gene_exprDF.loc[GENE_ID, bottom_df.columns.tolist() + top_df.columns.tolist()] 	# get expression values for only the top and bottom IDs

		# plot these expression levels as a bar plot
		plot_df = plot_data.to_frame()
		print(f'plot_df:\n{plot_df}')

		plot_df['GTEx samples'] = 'Top 20% OPA1'
		plot_df.loc[bottom_df.columns.tolist(),'GTEx samples'] = 'Bottom 20% OPA1'
		plot_df['GTEx sample ID'] = plot_df.index
		plot_df.sort_values(axis=0,by=GENE_ID,inplace=True)
	#	plot_df.to_csv(rf'GTEx_LV_OPA1_topBot20_29-07.csv')

	# 	quintile_expr = sns.barplot(data=plot_df, x='GTEx sample ID', y=GENE_ID, hue='GTEx samples', palette=COLOURS, dodge=False)
	# 	quintile_expr.get_xaxis().set_ticks([])
	# 	quintile_expr.set_xlabel('GTEx Left Ventricle - top & bottom OPA1 quintiles')
	# 	quintile_expr.set_ylabel('OPA1 expression (TMM-normalized TPM)')
	# 	quintile_expr.set_ylim(0.00,45)
	# 	# geneID_expr = expression_analysis.split_bar(plot_data=plot_df, plot_x='GTEx sample ID', plot_y=GENE_ID, plot_hue='GTEx samples')
	# 	# quintile_expr.text(0.0, 0.55, f'OPA1 LV mRNA (TPM, TMM normalized)', va='center', rotation='vertical') 	# Add y-axis text across both subplots
	# 	fig = quintile_expr.get_figure()
	# 	plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")	
	# 	fig.savefig(f'outputs/OPA_quintiles_LVexpr_15-08.png', bbox_inches='tight',dpi=600)
	# 	plt.close()
	# 	exit()
	# 	## As a sanity check, we would like to compare the overall gene expression distributions in our groups
	# 	## Compute significance tests on the distributions of expression values between each group (Mann-Whitney U test)
	# 	## Create histograms to understand the distribution of gene expression in our different groups

	# 	print(f'exp_df:\n{expr_df}')
	# 	print(f'top cols:\n{top_df.columns.tolist()}\nBottom cols:\n{bottom_df.columns.tolist()}')

	# 	# Compute statistics
	# 	group1 = top_df.values.flatten()
	# 	group2 = bottom_df.values.flatten()
	# 	u_stat, u_pvalue = mannwhitneyu(group1, group2)
	# 	print(f'Mann-Whitney U test: Statistic={u_stat}, p-value={u_pvalue}')
	# 	exit()
	# 	# expr_plt = expression_analysis.plot_exprDistribution(expr_df=expr_df, 
	# 	# 	group_ids=[top_df.columns.tolist(),  bottom_df.columns.tolist()], 
	# 	# 	group_names=['Top 20% OPA1', 'Bottom 20% OPA1'], 
	# 	# 	group_colours=COLOURS[::-1] 	# gives nicer overlay of colours
	# 	# 	)

	# 	# pval = f'{u_pvalue:.6f}'
	# 	# expr_plt.set_xlabel('mRNA Expression (TPM, TMM normalized)')
	# 	# expr_plt.set_title(f'GTEx LV gene expression in OPA1 quintiles (p-val: {pval})')

	# 	# # Add lines indicating the median of each distribution
	# 	# ymin, ymax = expr_plt.get_ylim()
	# 	# print(f'top median: {np.median(group1)}, bottom {np.median(group2)}')
	# 	# print(f'top mean: {np.mean(group1)}, bottom {np.mean(group2)}')
	# 	# # expr_plt.vlines(x=[np.median(group1),np.median(group2)], ymin=0, ymax=ymax, colors=COLOURS[1], ls='--', lw=1)

	# 	# fig = expr_plt.get_figure()
	# 	# fig.savefig(f'outputs/LV_heart_GTExexpr_hist_loglog_TMM_29-07.png', dpi=600)
	# 	# plt.close()

	# 	# We want to now compare pathologies of these samples
	# 	# Get tissue samples for each expression IDs
	# 	# Split pathology notes in to list and check for counts of different notes in our high/low groups and overal
	# 	# Get proportion (probability) for each group
	# 	# Can plot as bars or some kind of probability plots
	# # plot_df.to_csv(rf'new_high_low_expressors_OPA_vent_allAges.csv')

	# # A final sanity check is to compare expression between "old" and "young" samples we are using
	# # This assumes that all, or most age brackets, 20-29,...70-79, in GTEx were used (AGE_CUTOFF is large).
	# # We split ages a <50 and =>50, if a lower age cutoff it used this may be skewed.
	# young_ids = format_sampleDF.loc[format_sampleDF.loc[:,'Age'] < 50]
	# old_ids = format_sampleDF.loc[format_sampleDF.loc[:,'Age'] > 49]
	# print(f'young_ids:\n{young_ids["Tissue Sample ID"]}')
	# print(f'top_df:\n{top_df.columns}\n{young_ids["Tissue Sample ID"].to_list()}')

	# young_top = expression_analysis.subset_expr(subset_ids=young_ids, expr_df=top_df)
	# print(f'young_top:\n{young_top}')
	# old_top = expression_analysis.subset_expr(subset_ids=old_ids, expr_df=top_df)
	# print(f'old_top:\n{old_top}')

	# young_bottom = expression_analysis.subset_expr(subset_ids=young_ids, expr_df=bottom_df)
	# print(f'young_bottom:\n{young_bottom}')
	# old_bottom = expression_analysis.subset_expr(subset_ids=old_ids, expr_df=bottom_df)
	# print(f'old_bottom:\n{old_bottom}')

	# dfs = [young_bottom, old_bottom, young_top, old_top]
	# opa_dfs = [df.loc[GENE_ID] for df in dfs]

	# opa_df = pd.DataFrame(data=opa_dfs, index=['Young bottom 20%', 'Old bottom 20%', 'Young top 20%', 'Old top 20%'])
	# print(f'dfs:\n{dfs}\n\nopa_dfs:\n{opa_dfs}')
	# print(f'opa_df:\n{opa_df}')

	# # We want to distinguish the old and young subgroups, which there are several ways to do
	# # One option, is to differentiate them by colour
	# bottom_rgb = mcolors.ColorConverter.to_rgb(COLOURS[0])
	# top_rgb = mcolors.ColorConverter.to_rgb(COLOURS[1])

	# # Convert RGB to HLS, where we can adjust the lightness value
	# b_h, b_l, b_s = colorsys.rgb_to_hls(*bottom_rgb)
	# t_h, t_l, t_s = colorsys.rgb_to_hls(*top_rgb)

	# # Adjust lightness and convert back to RGB
	# l_bottom_rgb = colorsys.hls_to_rgb(b_h, b_l*1.5, b_s)
	# l_top_rgb = colorsys.hls_to_rgb(t_h, t_l*1.5, t_s)

	# # Create our palette with these colours
	# plt_palette = {
	# 	'Young top 20%': l_top_rgb, 
	# 	'Old top 20%': top_rgb, 
	# 	'Young bottom 20%': l_bottom_rgb, 
	# 	'Old bottom 20%': bottom_rgb
	# }
	# # Want to plot the distributions of expression in each group, and compare statistically

	# expr_plt = sns.violinplot(
	# data=opa_df.T,
	# palette=plt_palette,
	# saturation=0.85,
	# inner="point",
	# linewidth=1,
	# inner_kws=dict(linewidth=10)
	# )

	# expr_plt.set_title('OPA1 expression in quintiles split by age')
	# expr_plt.set_ylabel('OPA1 expression (TMM-normalized TPM)')
	# expr_plt.set_xlabel('OPA1 expression groups')
	# plt.xticks(rotation=30)

	# max_expr = opa_df.max().max()	# Get max expr for adjust ylims
	# plt.ylim(0, max_expr*1.2) 
	# fig = expr_plt.get_figure()
	# fig.savefig(f'outputs/LV_heart_GTExexpr_AgedOPA_comp_viol_31-07.png', bbox_inches='tight', dpi=600)
	# plt.close()

	# # Compute statistics between our aged subsets
	# top_t, top_pval = ttest_ind(opa_dfs[2], opa_dfs[3], nan_policy='omit')
	# print(f'T-test results for top groups: \n{opa_df}\n{opa_df}\n')
	# print(f'Statistics & p-value: {top_t}, {top_pval}')

	# bottom_t, bottom_pval = ttest_ind(opa_dfs[0], opa_dfs[1], nan_policy='omit')
	# print(f'T-test results for bottom groups: \n{opa_df}\n{opa_df}\n')
	# print(f'Statistics & p-value: {bottom_t}, {bottom_pval}')

	# for group in opa_dfs:
	# 	print(f'summary stats:\n{group.describe()}')
	# exit()
	# Subset our top and bottom DFs by age
	# young_top = top_df.loc[GENE_ID, young_ids["Tissue Sample ID"]]

	'''
	To further compare our groups, get differentially expressed genes between our groups
	Define differential expression as median log2FC > |1| and BH-adjusted Wilcoxon p-values < 0.05
	'''
	# Get differentially expressed genes between high/low groups
	groups_exprDF = gene_exprDF.loc[:, bottom_df.columns.tolist() + top_df.columns.tolist()] 	# get expression values for only the top and bottom IDs (NOTE: reverse bottom list to be ordered largest-> smallest)

	p_vals = groups_exprDF.apply(lambda x: expression_analysis.gene_rank_sum(geneExpr=x, group_a=top_df.columns.tolist(), group_b=bottom_df.columns.tolist()), axis=1)

	groups_exprDF['wilcoxon p-values'] = p_vals
	groups_exprDF['BH adjusted p-values'] = fdrcorrection(p_vals.values, method='indep', alpha=0.05)[1] 	# get adjust p-values

	sig_genes = groups_exprDF.where(groups_exprDF['BH adjusted p-values'] < 0.05)

	fc = expression_analysis.fold_change(geneExpr=groups_exprDF, group_a=top_df.columns.tolist(), group_b=bottom_df.columns.tolist(), method='median', log2=True)
	groups_exprDF.loc[:,'Log2 fold change'] = fc
	groups_exprDF.to_csv(rf'test_sigOPA_genes_take2.csv')
	exit()

	# Now get significantly differentially erxpressed genes with FC > |1.|5
	sigDEGs = groups_exprDF.where((groups_exprDF['BH adjusted p-values'] < 0.05) & ((groups_exprDF['Log2 fold change'] <= -1) | (groups_exprDF['Log2 fold change'] >= 1)))
	print(f'sigDEGs shape: {sigDEGs.shape}')
	sigDEGs.dropna(axis=0, how='all', inplace=True)
	print(f'sigDEGs non-zero shape: {sigDEGs.shape}')
	# sigDEGs.to_csv(rf'test_sigOPA_DEGs_1FCcsv')


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
	fig.savefig(f'outputs/OPA1_quintiles_volcano_29-07.png', bbox_inches='tight', dpi=600)
	plt.close()

	DE_plot = expression_analysis.volcano(data=groups_exprDF, 
			x='Log2 fold change', 
			y='BH adjusted p-values', 
			x_thresh=1.5, 	# specific our log2 FC threshold (absolute)
			y_thresh=0.01, 	# specify our FDR threshold
			sig_palette=COLOURS[::-1]
			)

	fig = DE_plot.get_figure()
	fig.savefig(f'outputs/OPA1_quintiles_volcano_inverse_29-07.png', bbox_inches='tight', dpi=600)
	plt.close()
	exit()
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
	fig.savefig('outputs/LV_OPA-PLSDA_29-07.png', dpi=400)
	exit()

	weights_df.to_csv(rf'outputs/test_plsda_weights_29-07-24.csv')
	exit()

	'''
	Similarily, we can use a Relief-based method for feature importance.
	Here, we can treat our data as a regression tast where we want to identify the most important
	genes for predicting OPA1 expression in our samples from top/bottom expressors
	'''

	target_expr = groups_exprDF.loc[GENE_ID].values 	# get OPA expression values for regression targets
	other_expr = groups_exprDF.drop(labels=GENE_ID, axis=0) 	# All other values as features
	print(f'other_expr shape: {other_expr.values.shape}\noother_expr[:-1] shape: {other_expr[:-1].values.shape}')
	gene_ids = other_expr.index.tolist()
	print(f'len(gene_ids): {len(gene_ids)}')
	# exit()

	relief_scores = expression_analysis.run_relief(data=other_expr.values.transpose(), 
		targets=target_expr, 
		discrete=False,
		multiround=True,
		multi_headers=gene_ids
		)

	for feature_name, feature_score in zip(gene_ids, relief_scores):
		print(f'{feature_name}:    {feature_score}')

	# # Cross-validation evaluation scores for non-multiround ReliefF (since TurF requires a lot of computation)
	# _, relief_evals = expression_analysis.run_relief(data=other_expr.values.transpose(), 
	# 	targets=target_expr, 
	# 	discrete=False,
	# 	multiround=False,
	# 	multi_headers=gene_ids,
	# 	evaluate=True
	# 	)

	# print(f'scores (train,test): {relief_evals}')
	# mean_evals = [np.mean(metrics) for metrics in relief_evals]
	# print(f'mean score: {mean_evals}')


	relief_df = pd.Series(data=relief_scores, index=gene_ids, name='Relief Scores')
	relief_df.sort_values(ascending=False, inplace=True)
	relief_df.to_csv(rf'relief_scores_k20_scaled.csv')
main()