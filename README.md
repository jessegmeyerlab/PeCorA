# Peptide Correlation Analysis (PeCorA)
is an analysis strategy to find peptides within proteins with quantities that differ significantly from other peptides in that protein across treatment groups. At a high level, PeCorA achieves this by looking for a difference in slope of each peptide quantity across treatment group versus all other peptides assigned to the same protein. Practically, this is achieved using a linear model to assess the statistical p-value of the interaction term of peptide group and treatment group. 

PeCorA is written in R and can be run from the command line as described below. 

### Preprint:
https://www.biorxiv.org/content/10.1101/2020.08.21.261818v2

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


### Output
Successful analysis will output the following to the directory specified by -d/--dir:
1. global_data_scaling.png, shows the distributions of peak areas per sample before and after scaling
2. allplots/, folder containing barplots that visualize the peptides with quantitative differences versus the other peptides in that protein
3. out.txt, table containing all the peptides that were tested, the proteins they belong to, the raw p-value, and the corrected p-value. 
