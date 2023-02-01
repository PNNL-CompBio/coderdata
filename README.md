## CANDLE Data Processing

This repository leverages existing datasets to collect the data required for CANDLE data analysis. Since each deep learning model requires distinct data capabilities, the goal of this repository is to collect and format all data into a schema that can be leveraged for existing models.

## Aggregated data tables

We will store the data in tables that are represented by the files below. Each data-specific model can be generated from a smaller set of these tables. The schema for these tables is represented below. 


The files will be comma-delimited and named follows:

1. genes.csv
2. drugs.csv
3. samples.csv
4. experiments.csv
5. expression.csv
6. mutations.csv
7. copynumber.csv


## PharmacoGX processing

One way to assemble these data is to use the PharmacoGX package and the curve fitting code. To do so, we have created custom scripts in the [pgx](pgx/) directory that collect the data available. Each dataset has slightly different data so we have collated it as needed and put it into the tables above.

### Gene mapping 
This set of scripts pulls data from pharmacoGX R package and formats each dataset to the schema above. 

### Sample mapping
The sample mapping file is derived from the DepMap identifiers and cellosaurus for now. 

### Curve fitting
The curve fitting is all run through a modified version of the code at the [curve repository](). 

## Model specific files

We then need to capture, for each model, the data necessary to train/build model. The code to do this should operate on the standard schema defined above so it can be re-run as we collect new data

Once we get a basic schema assigned we can start to write this code.
