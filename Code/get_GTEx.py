"""
Luke Kennedy - Fong-McMaster et al. (2024)

Code for subsetting bulk GTEx RNA-Seq data for analysis as reported in Fong-McMaster et al. (2024)
	Given a transcriptomics file for all GTEx samples (Download: https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression)
	we want to subset the data for our samples of interest (e.g., a tissue type) and save the subset in a smaller, more convenient file.
"""

# import modules
import pandas as pd
import numpy as np


"""
Given a list of GTEx sample IDs, we want to extract data for each sample and save them as outputs.
GTEx sample IDs correspond to, for example, samples with high or low muscle SLC7A11 mRNA. Want to call API to get all muscle MRNA for all genes in these individuals to identify other differentially expressed genes.
"""

# Define a method to get GTEx sample IDs for specified subejcts and/or tissues
# extract sample IDs for specific tissue samples for subjects from 'GTEx_histologyNotes.csv', which is used to extract gene expression data 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'
def get_sampleIDs(sample_info, subject_ids=None, tissues=None):

	# load file containing GTEX IDs, GTEx_histologyNotes.csv, which has subject IDs and the tissue sample IDs we need
	info_df = pd.read_csv(sample_info)

	# Check if any subject IDs are specified and extract rows with those IDs
	if isinstance(subject_ids, str):
		info_df = info_df.loc[info_df['Subject ID'] == subject_ids]

	elif isinstance(subject_ids, list):
		info_df = info_df.loc[info_df['Subject ID'].isin(subject_ids)]

	# Check if any tissue types are specified and extract rows for those tissues
	if isinstance(tissues, str):
		info_df = info_df.loc[info_df['Tissue'] == tissues]

	elif isinstance(tissues, list):
		info_df = info_df.loc[info_df['Tissue'].isin(tissues)]

	# Now have our desired subset, we want to extract just the Tissue Sample IDs for use in gene expression file
	return info_df['Tissue Sample ID'].tolist()


# Define a method to search the experssion data for a list of specified IDs/tissues
# GTEx gene expression data is very large (~55000 x 18000), usually we only want to look at gene expression in some tissues and/or subjects. This method will load only those data from the gene expression file
def get_subjectExpr(gtex_expr, gtex_samples, gtex_subjectIDs=None, gtex_tissues="Lung"):

	# Get the list of IDs to look for in our expression data
	if gtex_subjectIDs is None and gtex_tissues is None:
		print('\n\n!!!!\nNo data specified, loading ALL GTEX EXPRESSION DATA (4GB)\n!!!!!\n\n') 	# Warning in case we try loading all data

	select_ids = get_sampleIDs(subject_ids=gtex_subjectIDs, tissues=gtex_tissues, sample_info=gtex_samples) 	# get tissue sample IDs

	# GTEx expression data have additional ID information than what is in other datasets
	# Can look for our IDs as subsets of GTEX IDs (NOTE: possibly some IDs do not overlap)

	# Load expression data, but only the column names
	expr_IDs = pd.read_csv(gtex_expr, header=None, delimiter='\t', skiprows=2, nrows=1) 	# import the first (non-meta data) column, which contains all sample ID names

	expr_IDs = expr_IDs.iloc[0,:].tolist() 	# convert to list for iterating

	# Now check our gene expression IDs for our selected IDs
	select_exprCols = ['Name', 'Description'] 	# This list is the columns of expression data we want, so Gene name and info + all of our sample IDs of interest
	for sample_ID in select_ids:
		matched = [expr_id for expr_id in expr_IDs if sample_ID in expr_id]
		if len(matched) != 1:
			print(f'For tissue Sample ID: {sample_ID}, {matched} are matched in expression data')

		select_exprCols.extend(matched)
	

	# Now want to load our data of interest
	expr_df = pd.read_csv(gtex_expr, delimiter='\t', skiprows=2, usecols=select_exprCols) 	# import the first (non-meta data) column, which contains all sample ID names

	return expr_df