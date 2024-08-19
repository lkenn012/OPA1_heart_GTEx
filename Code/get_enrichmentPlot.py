"""
Code for generating enrichment plots of ShinyGO data for enrichments in genes correlated to Slc7a11 in C2C12 differentiation datasets as reported in Kanaan et al (2023)
	- Loads a data file contain ShinyGO ouputs for both positive and negative correlated gene sets
	- Transforms values (i.e. FDR values) for better visualization
	- Plots and formats the figures
	- Save figure

"""

# import modules
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


"""
load our enrichment data, modified from ShinyGO (such that negative correlations are denoted by negative fold enrichments)
	Biological process: corrs_enrichBioProc.csv
	Cellular Component: corrs_enrichCellCom.csv
	Molecular Function: corrs_enrichMolFunc.csv
	KEGG: corrs_enrichKEGG.csv
"""

# Define main function for 
def main():


	'''
	Arguments
	UP_FILE (str):			The path to a file containing a ShinyGO UP-regulated enrichment
	DOWN_FILE (str):		The path to a file containing a ShinyGO DOWN-regulated enrichment

	COLOURS (list, opt):	Extreme values for creating a colour map for enrichments (RGBs)
	MAX_PATHS (int, opt):	Defines the maximum number of pathways to show (in each group)
	'''


	UP_FILE = r'ShinyGo\OPA_conserved_UP_KEGG.csv'
	DOWN_FILE = r'ShinyGo\OPA_conserved_DOWN_KEGG.csv'

	COLOURS = [[249,165,27], [58,84,164]] 	# Orange and blue
	MAX_PATHS = 10


	# Get our dataframes
	up_enrich = pd.read_csv(UP_FILE)
	down_enrich = pd.read_csv(DOWN_FILE)

	# Select the most enriched pathways
	if MAX_PATHS:
		up_enrich = up_enrich.nlargest(MAX_PATHS, columns='Fold Enrichment')
		down_enrich = down_enrich.nlargest(MAX_PATHS, columns='Fold Enrichment')

	down_enrich['Fold Enrichment'] = -1*down_enrich['Fold Enrichment'] 	# Assign negative enrichment for down-regulated genes

	# combine and format enrichments for plotting
	enrich_df = pd.concat([up_enrich, down_enrich])
	enrich_df.sort_values('Fold Enrichment', inplace=True)
	enrich_df['-log10(FDR)'] = -np.log(enrich_df.loc[:,'Enrichment FDR']) * np.sign(enrich_df.loc[:,'Fold Enrichment'])

	# Pathway names include a prefix (Path:hsa0000, GO:00000) which can be removed for clarity
	# Define a regex pattern to match and just keep the name of the paths
	pat = r'\w+\:\w+\s'
	temp = enrich_df['Pathway'].str.replace(pat,'',regex=True)
	enrich_df["Pathway"] = temp

	"""
	Actual code related to plotting the data
		Define the colourmap from COLOURS for shading our nodes
		Define layout and design of the figure
		plot
		Add missing features
		save
	"""

	# Define palette
	if COLOURS:
		# convert hexocdes to HSL values for generating palette
		lower_rgb = [val/255 for val in COLOURS[0]]	# mcolors expects RGB to be in range [0,1], standard RGB is [0,255]
		upper_rgb = [val/255 for val in COLOURS[1]]
		lower_hsv = mcolors.rgb_to_hsv(lower_rgb)
		upper_hsv = mcolors.rgb_to_hsv(upper_rgb)
		print(f'upper & lower hsvs: {lower_hsv}, {upper_hsv}')
		# Seaborn expects hue from HSV in degrees, not float [0,1]
		neg_hue = lower_hsv[0] * 360
		pos_hue = upper_hsv[0] * 360

		# Create a diverging palette
		palette = sns.diverging_palette(neg_hue, pos_hue, s=90, l=40, sep=1, as_cmap=True)
	else:
		palette = 'rdBlu_r'

	# Define the style, font, etc.
	sns.set_style('whitegrid')
	sns.set_context("paper", rc={"font.size":8,"axes.titlesize":16})
	plt.rcParams['font.family'] = 'Arial'
	fig, ax = plt.subplots(figsize=(6, 6))

	sns.scatterplot(data=enrich_df, 
		 x='Fold Enrichment',
		 y='Pathway',
		 size='nGenes', 
		 sizes=(min(enrich_df.loc[:,'nGenes'])*7, max(enrich_df.loc[:,'nGenes'])*7), 	# Range of size to map {size} to values in our data, adjust scaling so larger values do not overlap.
		 hue='-log10(FDR)', 
		 palette=palette, 
		 edgecolor='black', 	# edge color of points
		 ax=ax
		 )

	plt.xlim(-max(abs(enrich_df.loc[:,'Fold Enrichment']))-2, max(abs(enrich_df.loc[:,'Fold Enrichment'])+2)) 	# specify x-axis range
	plt.title('KEGG Pathway Enrichment', fontsize=14) 	# NOTE: replace with title of enrichment

	# Now formatting stuff
	# By default, legend will display points for size (nGenes) and hue (significance), however we will be using a color bar for hue so only want to include size values
	h, l = ax.get_legend_handles_labels()
	leg = ax.legend(h[8:], l[8:], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0, title='Genes', frameon=False)

	ax.axvline(0, linewidth=2, color='gray') 	# Gray line at 0 for reference
	ax.xaxis.grid(False)
	plt.tick_params(axis='x', bottom=True, color='gray') 	# Add ticks at x-axis
	ax.invert_yaxis()

	"""
	For continuous FDR values, need to create a colorbar for the legend
	"""

	# Define limits of the color bar based on largest absolute values
	if enrich_df.loc[:,'-log10(FDR)'].abs().max() >= enrich_df.loc[:,'-log10(FDR)'].abs().min():
		color_lim = enrich_df.loc[:,'-log10(FDR)'].abs().max()
	else:
		color_lim = enrich_df.loc[:,'-log10(FDR)'].abs().min()

	# Create colorbar map
	norm = plt.Normalize(color_lim, -color_lim)
	sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
	sm.set_array([])

	# Create axis object for colorbar and format
	cax = fig.add_axes([ax.get_position().x1+0.055, ax.get_position().y0+0.22, 0.075, ax.get_position().height/3])
	cb = ax.figure.colorbar(sm, fraction=0.04, cax=cax) #, drawedges=True)
	cb.outline.set_color('black')
	cax.tick_params(size=0)
	# cax.set_color('black')
	cax.set_xlabel('-log10(FDR)', horizontalalignment='left', x=-0.45)
	cax.xaxis.set_label_position('top') 

	# plt.xticks(ax.get_xticklabels()) 	# Add x-axis ticks
	ax.spines['bottom'].set_color('gray') 	# Adjust bottom border to anchor the figure
	ax.spines['bottom'].set_linewidth(2)

	# Save figure
	plt.savefig(rf'outputs\enrichment_KEGG2_test.tiff', bbox_inches='tight', pad_inches=0.9, dpi=600)

# RUN MAIN()
main()