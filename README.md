# Peptide Correlation Analysis (PeCorA)
is an analysis strategy to find peptides within proteins with quantities that differ significantly from other peptides in that protein across treatment groups. At a high level, PeCorA achieves this by looking for a difference in slope of each peptide quantity across treatment group versus all other peptides assigned to the same protein. Practically, this is achieved using a linear model to assess the statistical p-value of the interaction term of peptide group and treatment group. 

PeCorA is written in R and can be run from the command line as described below. 

### Usage 
Rscript.exe PeCorA.R 

Required:
-f filename.csv 
File containing table in long format of peptides, their quantities, and the proteins they belong to. This file must contain the following columns (check spelling and letter case):
1. "Condition" - group labels of the conditions. Can be more than 2 but must be at least 2. 
2. "Peptide.Modified.Sequence" - peptide sequence including any modifications
3. "BioReplicate" - numbering for biological replicates
4. "Protein" - protein membership for each peptide

Options:

-o --out
Output file name

-d --dir
Directory to write output table and images
example, D:\output\directory

-c --cont
Name of the control group within the "Conditions" column

-a --area 
Name of column with peptide peak areas

-p --pval
Threshold to use for significant p-value when printing plots and printing summary results to console numbers, default =0.01

