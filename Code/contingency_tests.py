### statistical analyses of contingency tables

# import modules
from scipy.stats import fisher_exact, chi2_contingency


# Define the tables to test
# arrays = [[a,b]  	# [group1_cond1, group1_cond2], [group2_cond1, group2_cond2]
# 		   [c,d]]
arrays = [[[5,77],[3,79]],
	[[6,76],[4,78]],
	[[21,61], [28,54]],
	[[14,68],[4,78]],
	[[44,36],[56,26]],
	[[7,75],[1,81]],
	[[12,70],[16,66]],
	]
arr_names = ['top_bottom_hyper', 'top_bottom_infarc', 'top_bottom_ischem', 'top_bottom_slow', 'top_bottom_vent', 'top_bottom_inter', 'top_bottom_fastNat'] 	# We can also give names for our tests if we want

arrays2 = [[[28,54],[209,552]],
	[[16,66],[222,539]],
	[[12,70],[222,539]],
	[[44,36],[343,418]],
	[[56,26],[343,418]]
	]
names2 = ['bottom_overall_ischem', 'bottom_overall_fastNat', 'top_overall_fastNat', 'top_overall_vent', 'bottom_overall_vent']
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
		fish_res, fish_pval = fisher_exact(arr)
		print(f'Fisher exact res & p-value: {fish_res}, {fish_pval}')

		chi_res = chi2_contingency(arr)
		print(f'Chi2 rest & p-value: {chi_res[0]}, {chi_res[1]}')

	return

cont_tests(arrays, arr_names)

print(f'\n!!! Overall tests   !!!\n')
cont_tests(arrays2,names2)