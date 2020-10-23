# Peptide correlation analysis (PeCorA) <img src="man/figures/PECORA_hex.png" align="right" height="150"/>

`PeCorA` is short for "Peptide Correlation Analysis", which is an R package that enables detection of discordant peptide quantities in shotgun proteomics data. The package also contains a number of published proteomics datasets processed with different processing tools to demonstrate the workflow.

Please find the relevant preprint here: https://doi.org/10.1101/2020.08.21.261818

### Install

You can install `PeCorA` from github downloading the package by cloning
the repository.

`$ git clone https://github.com/jessegmeyerlab/PeCorA.git`

`$ R CMD INSTALL PeCorA-master`

Alternatively you can install `PeCorA` directly from R using devtools:

``` r
library(devtools)
install_github("jessegmeyerlab/PeCorA")
```

Or you can install`PeCorA` from CRAN by typing in R:
install.packages(“PeCorA”)

Once installed, load the package by writing in the console

``` r
library(PeCorA)
```

### Available datasets

Currently, there are three datasets available in `PeCorA`.

| Data              | Description                                                                                                                                              |
| :---------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------- |
| PeCorA\_noZ       | Primary mouse microglia dataset (PXD014466)                                                                                                              |
| input.dda.iprg.pg | BRF Proteome Informatics Research Group (iPRG) 2015 Study: Detection of Differentially Abundant Proteins in Label-Free Quantitative LC-MS/MS Experiments |
| Covid\_peptides   | Large-scale proteomic Analysis of COVID-19 Severity                                                                                                      |

### Loading data

Data available in the package is loaded into the `R` session using the
`load` function; for instance, to get the DDA iPRG data from [Choi et
al 2017](https://pubmed.ncbi.nlm.nih.gov/27990823/):

``` r
data("input.dda.iprg.pg")
```

To get more information about a dataset, see its manual page.

``` r
?input.dda.iprg.pg
```

### How to use

PeCorA requires a filename.csv file containing table in long format of
peptides, their quantities, and the proteins they belong to. This file
must at least contain the following columns (check spelling and letter
case):

“Condition” - group labels of the conditions. Can be more than 2 but
must be at least 2. “Peptide.Modified.Sequence” - peptide sequence
including any modifications “BioReplicate” - numbering for biological
replicates “Protein” - protein membership for each peptide

You may need to transform your data into PeCorA-ready format. For
example ransform peptides.txt output of MaxQuant into t use function
`import_preprocessed_for_PeCorA`.

### Functions

The main function of the package is called `PeCorA`, which fits a linear
model with interaction between peptides and biological treatment groups.

### Vignette

See
[Vignette](https://github.com/jessegmeyerlab/PeCorA/blob/master/vignettes/PeCorA_vignette.pdf)
for complete workflow.

### Contact

If you have any questions or suggestions please contact us:

Maria Dermit : maria.dermit at qmul.ac.uk

Jesse Meyer: jesmeyer at mcw.edu
