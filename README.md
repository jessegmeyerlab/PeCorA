# Peptide Correlation Analysis (PeCorA)
is an analysis strategy to find peptides within proteins with quantities that differ significantly from other peptides in that protein across treatment groups. At a high level, PeCorA achieves this by looking for a difference in slope of each peptide quantity across treatment group versus all other peptides assigned to the same protein. Practically, this is achieved using a linear model to assess the statistical p-value of the interaction term of peptide group and treatment group. 

PeCorA is written in R and can be run from the command line as described below. 

### Preprint:
https://www.biorxiv.org/content/10.1101/2020.08.21.261818v2
To reproduce preprint results, use peptide quantities `data("PeCorA_noZ")` and see [Vignette](https://github.com/demar01/PeCorA/blob/master/vignettes/PeCorA_vignette.pdf) for complete workflow.

### Installation and usage 

Once installed, load the package by writing in the console

```{r}
library(PeCorA)
```
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

### Functions
The main function of the package is called `PeCorA`, which fits a linear model with interaction between peptides and biological treatment groups.

### Vignette 
See [Vignette](https://github.com/demar01/PeCorA/blob/master/vignettes/PeCorA_vignette.pdf) for complete workflow.

### Contact
If you have any questions or suggestions please contact us:

Maria Dermit : maria.dermit at qmul.ac.uk 

Jesse Meyer: jesmeyer at mcw.edu

