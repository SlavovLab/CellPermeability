# Limiting the impact of protein leakage in single-cell proteomics

This repository outlines the work completed as part of the study "Limiting the impact of protein leakage in single-cell proteomics" (Leduc et al. 2024).

## Tool for identifying permeablized cells in single-cell proteomics studies

If performing a single cell proteomics study, the use permeability stain such as Sytox Green can be helpful for excluding damage cells from downstream analysis. If permeablized cells were not excluded from the sample preparation process, they can still be identified computationally to avoid artifacts in your analysis. To facilitate this task, we developed a classifyer based on proteomic profiles from experimental data from intact and permeablized cells.

To use the classifyer, install the `QuantQC` package:

``` 
devtools::install_github("https://github.com/Andrew-Leduc/QuantQC")
library(QuantQC)
```

Once the package is installed, the relevant function is `FindPermeableCells()`. The input is a matrix (Protein X single cells) with Uniprot identifiers on the rownames. An optional argument is the species. The default is `species = "Mouse"`. Human is also supported, `species = "Human"`.

The output is a vector length number single cells. Each value is the probability of permeablization. The distribution may differ depending on the data set but generally values observed probabilities above 0.2 were observed as permeable cells. 


## Reproducing analysis from the paper

In the code directory there are several scripts. 

### Preprocess

QuantQC functions to do initial processing and QC of the LC-MS/MS data and metadata mapping. Outputs are Protein X single cell matricies for the Frozen and Fresh experiments and associated meta data files.

### Cell type asignment 

Protein X single cell matricies are integrated via the LIGER algorithm (Welch 2019) with mRNA-seq data from the same tissue to transfer cell type identification. Outputs are meta data files with an additional column for cell type identity.

### Downstream analysis

The final script takes in the meta data and protein data matricies and performs the analysis of protein abundance differences between permeable and intact cells from figures 1d,e and figure 2 from (Leduc et al 2024). 


### Cell type asignment 


## Link to study and additional resources

* **Preprint:** Leduc A, Xu Y, Shipkovenska G, Dou Z, Slavov N, [Limiting the impact of protein leakage in single-cell proteomics](https://www.biorxiv.org/content/10.1101/2024.07.26.605378v1), *bioRxiv*, doi: [10.1101/2024.07.26.605378](https://doi.org/10.1101/2024.07.26.605378)

* [Data Website](https://scp.slavovlab.net/Leduc_et_al_2024)

* [nPOP sample preparation Website](https://scp.slavovlab.net/nPOP) &nbsp; | &nbsp; [Download data from Leduc et al., 2022](https://scp.slavovlab.net/Leduc_et_al_2022) | &nbsp; [Download data from protocol](https://scp.slavovlab.net/Leduc_et_al_2023)

* [*Genome Biology* nPOP Article](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02817-5)
