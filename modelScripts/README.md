# Model scripts
This directory contains the scripts required to format the data from the harmonized dataset to the specific models at hand.

## Data Included
To date we have the datasets included in the [data](../data/) directory.

## Models
We are in the process of developing customized scripts for the following models.

| Model Name | Data Types | Script | External data/tools?|
| --- | --- | ---- | --- |
| [UNO]() | gene expression, mutation, cnv, proteomics, AUC/IC50, Drug fingerprints | | None}
| [DeepTTC](https://github.com/jianglikun/DeepTTC) | gene expression, AUC, SMILES | [./deep_ttc.py]()| None|
| [PathDSP](https://github.com/TangYiChing/PathDSP) | gene expression, AUC, mutation data | | GSEA analysis of expression, network analysis of mutations |
| [Paccman_MCA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7319576/) |IC50, drug smiles | | |
| [HiDRA](https://github.com/GIST-CSBL/HiDRA) |gene expression, morgan fingerprint of drug, IC50| TBD| None|
| [DrugGCN](https://www.mdpi.com/2227-7390/9/7/772) |gene expression, AUC, IC50|| protein interaction network, L1000?|
| [GraphDRP](https://github.com/hauldhut/GraphDRP)|gene expression, SMILES, |||
| [tCNNS](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2910-6) |IC50, gene expression, SMILES |||

## Additional tools
Some of these tools require additional helper tools such as GSEA or networks to implement. We will encode those functions here as well.
