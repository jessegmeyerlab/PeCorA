# Peptide Correlation Analysis (PeCorA)
is an analysis strategy to find peptides within proteins with quantities that differ significantly from other peptides in that protein across treatment groups. At a high level, PeCorA achieves this by looking for a difference in slope of each peptide quantity across treatment group versus all other peptides assigned to the same protein. Practically, this is achieved using a linear model to assess the statistical p-value of the interaction term of peptide group and treatment group. 

### Preprint:
https://www.biorxiv.org/content/10.1101/2020.08.21.261818v2

### Install
You can install `PeCorA` from github downloading the package by cloning the repository.

`$ git clone https://github.com/jessegmeyerlab/PeCorA.git`

`$ R CMD INSTALL PeCorA-master`

Alternatively you can install `PeCorA` directly from R using devtools:

```{r}
library(devtools)
install_github("jessegmeyerlab/PeCorA")
```

Or you can install`PeCorA` from CRAN by typing in R: `install.packages("PeCorA")`

Once installed, load the package by writing in the console:

```{r}
library(PeCorA)
```
To reproduce preprint results, use peptide quantities `data("PeCorA_noZ")` and 
see [Vignette](https://github.com/demar01/PeCorA-1/blob/master/vignettes/PeCorA_vignette.pdf) for complete workflow.

### Available datasets

Currently, there are three datasets available in `PeCorA`.

Data                    |Description                                                                                                          |
|:-----------------------|:--------------------------------------------------------------------------------------------------------------------|
|PeCorA_noZ   |Primary mouse microglia dataset (PXD014466)|
|input.dda.iprg.pg  |BRF Proteome Informatics Research Group (iPRG) 2015 Study: Detection of Differentially Abundant Proteins in Label-Free Quantitative LC-MS/MS Experiments |
|Covid_peptides  |Large-scale proteomic Analysis of COVID-19 Severity|

### Loading data

Data is loaded into the `R` session using the `load` function; for
instance, to get the DDA iPRG data from
[Choi et al 2017](https://pubmed.ncbi.nlm.nih.gov/27990823/):

```{r}
data("input.dda.iprg.pg")
```
### How to use
PeCorA requires a filename.csv file containing table in long format of peptides, their quantities, and the proteins they belong to. This file must at least contain the following columns (check spelling and letter case):

"Condition" - group labels of the conditions. Can be more than 2 but must be at least 2.\

"Peptide.Modified.Sequence" - peptide sequence including any modifications.\

"BioReplicate" - numbering for biological replicates.\

"Protein" - protein membership for each peptide.\

You may need to transform your data into PeCorA-ready format. For example ransform peptides.txt output of MaxQuant into t use function `import_preprocessed_for_PeCorA`.

### Functions
The main function of the package is called `PeCorA`, which fits a linear model with interaction between peptides and biological treatment groups.

### Vignette 
See [Vignette](https://github.com/demar01/PeCorA-1/blob/master/vignettes/PeCorA_vignette.pdf) for additional information.

### Contact
If you have any questions or suggestions please contact us:

Maria Dermit : maria.dermit at qmul.ac.uk 

Jesse Meyer: jesmeyer at mcw.edu

