# Information
Code and data for the analysis and plotting of OPA1 GTEx left ventricle transcriptomics transcriptomics data as reported in:  
Fong-McMaster, C., Pulente, S. M., Kennedy, L., Smith, T. K., Myers, S., Kanaan, M., Karam, C., Cope, M., Lorenzen-Schmidt, I., Goergen, C. J., Fullerton, M. D., Cuperlovic-Culf, M., Mulvihill, E. E., & Harper, M. E. (2024). OPA1 mediates cardiac function and metabolism: in silico and in vivo evidence. *bioRxiv*, 2024-08.  
❤️📄 [Read the paper here](https://doi.org/10.1101/2024.08.23.605375)

## Code
Contains all python code for the analysis of GTEx transcriptomics data and generation of imaages used in publication main text and supplementary figures.

## Data
GTEx data used in analysis and results files from ShinyGO used in generating enrichment plots

## Environment information
The following version of Python along with all accompanying libraries are used in the code provided here. All code was implemented on Windows with an installation time of approximately 15 minutes. 

Python   --3.8.18

### Scientific programming & machine learning packages
- numpy   --1.24.3
- pandas   --2.0.3
- scikit-learn   --1.3.0
- scipy   --1.10.1
- statsmodels   --0.14.0
- conorm	--1.2.0
- skrebate	--0.62

### Plotting packages
- matplotlib   --3.7.2
- seaborn   --0.12.2

Installation information:
```
conda install numpy
conda install pandas 
conda install scikit-learn   
conda install scipy
conda install statsmodels  
conda install matplotlib  
conda install seaborn  
pip install conorm
pip install skrebate

conda/pip install foo=1.2.3 		<- For installing specific versions
```

