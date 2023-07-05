## CANDLE Data Processing

This repository leverages existing datasets to collect the data
required for CANDLE data analysis. Since each deep learning model
requires distinct data capabilities, the goal of this repository is to
collect and format all data into a schema that can be leveraged for
existing models.

The goal of this repository is two-fold: First, it aims to collate and
standardize the data for all CANDLE related models. This requires
running a series of scripts to build and append to a standardized data
model. Second, it has a series of scripts that pull from the data
model to create model-specific data files that can be run by the data
infrastructure. 

# IMPROVE Data Model

The goal of the data model is to collate drug response data together with molecular data in a way that can be easily ingested by machine learning models. The overall schema is shown below.

<img src="origDataSchema.jpg" width=25% height=25%>

We will store the data in tables that are represented by the files below. Each data-specific model can be generated from a smaller set of these tables. The schema for these tables is represented below. 

The files are comma-delimited and named follows:
1. genes.csv
2. drugs.tsv.gz --> Drug names have commas and quotes in them, therefore require tab delimited
3. samples.csv
4. experiments.csv.gz --> compressed to fit on github
5. transcriptomics.csv.gz
6. mutations.csv.gz 
7. copy_number.csv.gz
8. methylation.csv.gz
9. mirnas.csv.gz

## Building the data model
The data model requires four steps, that need to be run in order
because they depend on each other.

| Data model step | Description/Dependencies| Script |
| --- | --- | --- |
| Build gene table | Uses the BioConductor Entrez Gene DB Annotations to map identifiers| [data/initialGeneDB.R](./data/initialGeneDB.R) |
| Build sample table | Uses DepMap and Cellosaurus to initialize cell line information. Other samples added later. | [data/initialSampleDB.R](./data/initialSampleDB.R) |
| Add in drug data, gene-based data | drug data will run through PharmacoGX, being added on an as-needed basis | [processPGXData.py](./processPGXData.py)
| --- | --- | --- |

### PharmacoGX processing

One way to assemble these data is to use the PharmacoGX package and the curve fitting code. To do so, we have created custom scripts in the [pgx](pgx/) directory that collect the data available. Each dataset has slightly different data so we have collated it as needed and put it into the tables above.

#### Gene mapping 
This set of scripts pulls data from pharmacoGX R package and formats each dataset to the schema above. 

#### Sample mapping
The sample mapping file is derived from the DepMap identifiers and cellosaurus for now. 

#### Curve fitting
The curve fitting is all run through a modified version of the code at the [curve repository](https://github.com/levinas/curve). 

# Model-specific data files

We then need to capture, for each deep learning model, the data necessary to train/build model. The code to do this should operate on the standard schema defined above so it can be re-run as we collect new data

Once we get a basic schema assigned we can start to write this code.
