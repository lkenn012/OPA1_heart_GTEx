"""
Luke Kennedy - Fong-McMaster et al. (2024)

Statistical tests applied to pathologies and types of death found in GTEx_histologyNotes.csv for our samples of interest.
	Arrays contain contingency tables, and lists of names describe the respective contingency tables.
	These data were obtained manually based on sample IDs determined in `get_expression_groups.py` and
	the pathologies in the histology notes files.

	Note: other pathologies are present, and were tested, though no significant differences were identified

"""

# import modules
from scipy.stats import chi2_contingency


# define method for computing and returning the results of chi2 contingencincy tests
def cont_tests(arrays, arr_names):
	for i, arr in enumerate(arrays):
		for a in arr:
			print(f'array: {a}, sum ={sum(a)}')
		if arr_names != []:
			try:
				print(f'Test name: {arr_names[i]}')
			except:
				print(f'No name')
		print(f'For arr = {arr}, stat tests:')
		chi_res = chi2_contingency(arr)
		print(f'Chi2 rest & p-value: {chi_res[0]}, {chi_res[1]}')

	return

# define main method to run chi2 tests on our input data
def main():

	'''
	Arguments
	GROUP_COMPARISONS  (list, strs): 			Contains descriptions for each contingency table in 'GROUP_CONTINGENCIES' for top group and bottom group comparisons, pathology category listed.
	GROUP_CONTINGENCIES  (list, ints): 			Contains pathology occurences for groups being compared for top and bottom group comparisons
	
	OVERALL_COMPARISONS  (list, strs): 			Contains pathology occurences for groups being compared for top/bottom group and overall sample comparisons
	OVERALL_CONTINGENCIES  (list, ints): 		Contains descriptions for each contingency table in 'GROUP_CONTINGENCIES' top/bottom group and overall sample comparisons, pathology category listed.
	''' 

	GROUP_COMPARISONS = ['top_bottom_hypertensions', 'top_bottom_infarction', 'top_bottom_ischemia', 
		'top_bottom_slowDeath', 'top_bottom_ventilator', 'top_bottom_intermediateDeath', 'top_bottom_fastNaturalDeath']
	GROUP_CONTINGENCIES = [
		[[5,77],[3,79]],
		[[6,76],[4,78]],
		[[21,61],[28,54]],
		[[14,68],[4,78]],
		[[44,36],[56,26]],
		[[7,75],[1,81]],
		[[12,70],[16,66]],
		]

	OVERALL_COMPARISONS = ['bottom_overall_ischemia', 'bottom_overall_fastNaturalDeath', 
		'top_overall_fastNaturalDeath', 'top_overall_ventilator', 'bottom_overall_ventilator']
	OVERALL_CONTINGENCIES = [
		[[28,54],[209,552]],
		[[16,66],[222,539]],
		[[12,70],[222,539]],
		[[44,36],[343,418]],
		[[56,26],[343,418]]
		]

	# Run tests on our top vs. bottom comparisons
	print(f'Results for comparisons of Top OPA1 expression group vs. Bottom expression group:\n')
	cont_tests(GROUP_CONTINGENCIES, GROUP_COMPARISONS)

	# Run tests for comparisons of our groups to overall sample pathology occurences
	print(f'\nResults for comparisons of OPA1 expression groups vs. overall samples:\n')
	cont_tests(OVERALL_CONTINGENCIES, OVERALL_COMPARISONS)

	return


# RUN MAIN
main()