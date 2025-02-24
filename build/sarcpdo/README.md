# Sarcoma Organoid Dataset

This dataset comes from a [recent publication](https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(24)00296-0?uuid=uuid%3Abb17ce2a-7cc0-48bd-a56d-3bc46a4d5541) with over 100 organoids
from patient sarcomas. 

## Samples

The sample information is contained in [this synapse table](https://www.synapse.org/Synapse:syn61894699/tables/)


## Omics data

Omics data is spread across two Synapse instances. Genomic mutation data can be found in [this synapse table](https://www.synapse.org/Synapse:syn61894695/tables/) and RNAseq data can be found in [this synapse table](https://www.synapse.org/Synapse:syn64333318/tables/).

## Drug information

Drug information can be found at [this synapse link](https://www.synapse.org/Synapse:syn61892224/tables/).

## Experimental measurements

Full experimental data (including complete AUC information) was not available for this multi-omics dataset, so we used the available viability scores in the drug information listed above (syn61892224), and labelled them as 'published auc' in the final experiments table. 
