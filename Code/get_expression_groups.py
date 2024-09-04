"""
Luke Kennedy - Fong-McMaster et al. (2024)

Code for accessing and analyzing bulk GTEx RNAseq data for evaluation of OPA1 expression as reported in Fong-McMaster et al. (2024)
	Given a gene expression data file (bulk or tissue-specific), this code selects top and bottom expressing samples for a gene of interest.
	Certain selection criteria can be applied through the argument variables ("AGE_CUTOFF", etc.) and 
"""


# import modules

import get_GTEx  	# helper functions
import expression_analysis 		# Helper functions

import numpy as np
import pandas as pd
import conorm 		# TMM-normalization
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

	BULK_GTEx = r'data\GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'  	# !!PLACEHOLDER!! This data file can be downloaded via: 
																					# 		https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
	GTEx_INFO = r'data\GTEx_histologyNotes.csv'  	# !!PLACEHOLDER!! Replace 'data' with your directory
	GTEx_TISSUE = 'Heart - Left Ventricle'

	GENE_ID = 'ENSG00000198836.8' 	# OPA1
	AGE_CUTOFF = 79
	CND_CONTROL = {'Hardy Scale': 'Ventilator case'}
	SEX_CONTROL = True
	N_CUTOFF = 24 		# 24 = OPA1 10% (x2) after outliers, vents only, all ages
	COLOURS = ['#F9A51B', '#3A54A4'] 	# Orange and blue

	# Get our data of interest, left ventricle, from bulk GTEx expression data
	# NOTE: Rather than loading the bulk expression data everytime, it may be easier to save the tissue RNAseq data and load that file for future use.
	temp = get_GTEx.get_subjectExpr(gtex_expr=BULK_GTEx, gtex_samples=GTEx_INFO, gtex_tissues=GTEx_TISSUE)

	# Save the file for further downstream analyses
	print(f'temp df for just left ventricle data:\n{temp}\nsahpe:{temp.shape}')
	temp.to_csv(rf'leftVent_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.csv')


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
			cnd_ids = format_sampleDF.loc[format_sampleDF.loc[:, key] == values]
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

	elif not SEX_CONTROL:
		top_ids, bottom_ids = expression_analysis.get_expressors(expr_df=gene_exprDF, n_samples=N_CUTOFF, sort_by=GENE_ID, clean_outliers=False)

		top_df = gene_exprDF.loc[:,top_ids]
		bottom_df = gene_exprDF.loc[:,bottom_ids]

	print(f'summary satistics for combined groups:\nTop:\n{top_df.describe()}')
	print(f'Bottom:\n{bottom_df.describe()}')

	# Plot bar plots of our high/low expressors
	plot_data = gene_exprDF.loc[GENE_ID, bottom_df.columns.tolist() + top_df.columns.tolist()] 	# get expression values for only the top and bottom IDs

	# plot these expression levels as a bar plot
	plot_df = plot_data.to_frame()

	plot_df['GTEx samples'] = 'Top 20% OPA1'
	plot_df.loc[bottom_df.columns.tolist(),'GTEx samples'] = 'Bottom 20% OPA1'
	plot_df['GTEx sample ID'] = plot_df.index
	plot_df.sort_values(axis=0,by=GENE_ID,inplace=True)

	##
	## Before proceeding, we want to now save the top and bottom expression files for use
	## in differential gene expression code - `diff_gene_expression.py`
	##
	top_df.to_csv(rf'top_LV_OPA1_RNAseq.csv')
	bottom_df.to_csv(rf'bottom_LV_OPA1_RNAseq.csv')
	plot_df.to_csv(rf'GTEx_LV_OPA1_topBot20.csv') 	# Save ID info

	quintile_expr = sns.barplot(data=plot_df, x='GTEx sample ID', y=GENE_ID, hue='GTEx samples', palette=COLOURS, dodge=False)
	quintile_expr.get_xaxis().set_ticks([])
	quintile_expr.set_xlabel('GTEx Left Ventricle - top & bottom OPA1 quintiles')
	quintile_expr.set_ylabel('OPA1 expression (TMM-normalized TPM)')
	quintile_expr.set_ylim(0.00,45)

	# # NOTE: If expression levels are dramatically different, we can use the `split_bar()` method to clarify the plots
	# geneID_expr = expression_analysis.split_bar(plot_data=plot_df, plot_x='GTEx sample ID', plot_y=GENE_ID, plot_hue='GTEx samples')
	# quintile_expr.text(0.0, 0.55, f'OPA1 LV mRNA (TPM, TMM normalized)', va='center', rotation='vertical') 	# Add y-axis text across both subplots

	fig = quintile_expr.get_figure()
	plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")	
	fig.savefig(f'OPA1_LVquintiles_bar.png', bbox_inches='tight',dpi=600)
	plt.close()

	## As a sanity check, we would like to compare the overall gene expression distributions in our groups
	## Compute significance tests on the distributions of expression values between each group (Mann-Whitney U test)
	## Create histograms to understand the distribution of gene expression in our different groups

	# Compute statistics
	group1 = top_df.values.flatten()
	group2 = bottom_df.values.flatten()
	u_stat, u_pvalue = mannwhitneyu(group1, group2)
	print(f'Mann-Whitney U test: Statistic={u_stat}, p-value={u_pvalue}')

	expr_plt = expression_analysis.plot_exprDistribution(expr_df=expr_df, 
		group_ids=[top_df.columns.tolist(),  bottom_df.columns.tolist()], 
		group_names=['Top 20% OPA1', 'Bottom 20% OPA1'], 
		group_colours=COLOURS[::-1] 	# gives nicer overlay of colours
		)

	pval = f'{u_pvalue:.6f}'
	expr_plt.set_xlabel('mRNA Expression (TPM, TMM normalized)')
	expr_plt.set_title(f'GTEx LV gene expression in OPA1 quintiles (p-val: {pval})')

	fig = expr_plt.get_figure()
	fig.savefig(f'LV_heart_GTExexpr_hist_loglog_TMM.png', dpi=600)
	plt.close()


	# A final sanity check is to compare expression between "old" and "young" samples we are using
	# This assumes that all, or most age brackets, 20-29,...70-79, in GTEx were used (AGE_CUTOFF is large).
	# We split ages a <50 and =>50, if a lower age cutoff it used this may be skewed.
	young_ids = format_sampleDF.loc[format_sampleDF.loc[:,'Age'] < 50]
	old_ids = format_sampleDF.loc[format_sampleDF.loc[:,'Age'] > 49]

	young_top = expression_analysis.subset_expr(subset_ids=young_ids, expr_df=top_df)
	old_top = expression_analysis.subset_expr(subset_ids=old_ids, expr_df=top_df)

	young_bottom = expression_analysis.subset_expr(subset_ids=young_ids, expr_df=bottom_df)
	old_bottom = expression_analysis.subset_expr(subset_ids=old_ids, expr_df=bottom_df)

	dfs = [young_bottom, old_bottom, young_top, old_top]
	opa_dfs = [df.loc[GENE_ID] for df in dfs]

	opa_df = pd.DataFrame(data=opa_dfs, index=['Young bottom 20%', 'Old bottom 20%', 'Young top 20%', 'Old top 20%'])

	# We want to distinguish the old and young subgroups, which there are several ways to do
	# One option, is to differentiate them by colour
	bottom_rgb = mcolors.ColorConverter.to_rgb(COLOURS[0])
	top_rgb = mcolors.ColorConverter.to_rgb(COLOURS[1])

	# Convert RGB to HLS, where we can adjust the lightness value
	b_h, b_l, b_s = colorsys.rgb_to_hls(*bottom_rgb)
	t_h, t_l, t_s = colorsys.rgb_to_hls(*top_rgb)

	# Adjust lightness and convert back to RGB
	l_bottom_rgb = colorsys.hls_to_rgb(b_h, b_l*1.5, b_s)
	l_top_rgb = colorsys.hls_to_rgb(t_h, t_l*1.5, t_s)

	# Create our palette with these colours
	plt_palette = {
		'Young top 20%': l_top_rgb, 
		'Old top 20%': top_rgb, 
		'Young bottom 20%': l_bottom_rgb, 
		'Old bottom 20%': bottom_rgb
	}

	# Want to plot the distributions of expression in each group, and compare statistically
	expr_plt = sns.violinplot(
	data=opa_df.T,
	palette=plt_palette,
	saturation=0.85,
	inner="point",
	linewidth=1,
	inner_kws=dict(linewidth=10)
	)

	expr_plt.set_title('OPA1 expression in quintiles split by age')
	expr_plt.set_ylabel('OPA1 expression (TMM-normalized TPM)')
	expr_plt.set_xlabel('OPA1 expression groups')
	plt.xticks(rotation=30)

	max_expr = opa_df.max().max()	# Get max expr for adjust ylims
	plt.ylim(0, max_expr*1.2) 
	fig = expr_plt.get_figure()
	fig.savefig(f'LV_OPA_AgeGroups_violin.png', bbox_inches='tight', dpi=600)
	plt.close()

	# Compute statistics between our aged subsets
	top_t, top_pval = ttest_ind(opa_dfs[2], opa_dfs[3], nan_policy='omit')
	print(f'T-test results for top groups: \n{opa_df}\n{opa_df}\n')
	print(f'Statistics & p-value: {top_t}, {top_pval}')

	bottom_t, bottom_pval = ttest_ind(opa_dfs[0], opa_dfs[1], nan_policy='omit')
	print(f'T-test results for bottom groups: \n{opa_df}\n{opa_df}\n')
	print(f'Statistics & p-value: {bottom_t}, {bottom_pval}')

	for group in opa_dfs:
		print(f'summary stats:\n{group.describe()}')

main()