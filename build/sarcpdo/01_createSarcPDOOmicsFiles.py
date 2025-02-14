import synapseclient
#import pyarrow
import pandas as pd
import numpy as np
#import polars as pl
import argparse
import os




def download_and_format_transcriptomic(synLoginObject, genesTable, samplesTable):

    transcriptomic_query = synLoginObject.tableQuery("select * from syn64333318")
    transcriptomic_data = transcriptomic_query.asDataFrame()
    # remve samples not meant to be included, if present
    samples_to_remove = ['SARC0128_Tumor', 'SARC0128_Organoids', 'SARC0129_Tumor', 'SARC0129_Organoids', 'SARC0130_Tumor', 'SARC0130_Organoids']
    matching_columns = transcriptomic_data.filter(items=samples_to_remove)
    transcriptomic_data.drop(matching_columns, axis=1, inplace=True)
    mapped_df = transcriptomic_data.merge(genesTable, left_on='gene_name', right_on='gene_symbol', how='left')

    # get names of columns with transcriptomic data
    value_variables = transcriptomic_data.columns[transcriptomic_data.columns.str.contains("SARC")]
    melted = mapped_df.melt(id_vars='entrez_id', value_vars = value_variables, var_name='sample_name')
    # drop null entrez_ids
    meltedNoNulls = melted[melted['entrez_id'].notna()]
    # drop duplicated rows
    #meltedNoNulls = meltedNoNulls.drop_duplicates()
    # fixing names that are '..Organoids' (plural) in transcriptomic column names but '..Organoid' (single) in sample sheet
    meltedNoNulls['edited_sample_name'] = meltedNoNulls['sample_name'].str.split("s", expand=True)[0]
     # join with 'improve_sample_id' in sample sheet
    samplesTable['other_id_no_dash'] = samplesTable['other_id'].str.replace("-2", "_2")
    melted_joined = meltedNoNulls.merge(samplesTable[['other_id', "other_id_no_dash", "improve_sample_id"]], 
                                        left_on='edited_sample_name', right_on='other_id_no_dash', how='left')
    melted_joined.drop('other_id', axis=1, inplace=True)
    melted_joined['source'] = "Synapse"
    melted_joined['study'] = "Landscape of Sarcoma"
    # select desired columns  - 'entrez_id', 'improve_sample_id', 'transcriptomics' 'source', study
    melted_joined_renamed =melted_joined.rename({'value' : 'transcriptomics', 'edited_sample_name' : 'other_id'}, axis=1)
    final = melted_joined_renamed[['entrez_id', 'improve_sample_id', 'transcriptomics', 'source', 'study']]
    #dropduplicates (see a few lines above - should be down here)
    final = final.drop_duplicates()
    return final

def download_and_format_genomic_mutation(synLoginObject, genesTable, samplesTable):
    mutationQuery=synLoginObject.tableQuery("select * from syn61894695")
    mutationDF = mutationQuery.asDataFrame()
    mutationDF['Sample_ID_Tumor'] = mutationDF['Sample_ID'] + "_Tumor"
    # left join with genes table 
    mutation_merged = mutationDF.merge(genes, left_on='Gene', right_on='gene_symbol', how='left')
   # drop null entrez_ids
    mutation_merged[mutation_merged['entrez_id'].isna()]
    #split gene name to include portion without exon
    mutation_merged["Name"] = mutation_merged["Name"].str.split("[ \(|]", expand=True)[0]
    # reformat variant classification column to be accepted by linkML and correct
    mutation_merged["variant_classification"] =mutation_merged['Canonical_Variant_Classification']

    #mutation_merged['variant_classification'] = 
    #mutation_merged['variant_classification'].replace("Missense", "Missense_Mutation", inplace=True)
    mutation_merged.replace({'variant_classification': "Missense"}, "Missense_Mutation", inplace=True)
    #mutation_merged['variant_classification'] = 
    mutation_merged.replace({'variant_classification': "Splice_Donor"}, "Splice_Site", inplace=True)
    mutation_merged.replace({'variant_classification': "Splice_Acceptor"}, "Splice_Site", inplace=True)
    mutation_merged.replace({'variant_classification': "Nonsense"}, "Nonsense_Mutation", inplace=True)
    mutation_merged.replace({'variant_classification': "intron"}, "Intron", inplace=True)
    mutation_merged.replace({'variant_classification': "synonymous"}, "Silent", inplace=True)
    mutation_merged.replace({'variant_classification': "Inframe_Del"}, "In_Frame_Del", inplace=True)
    mutation_merged.replace({'variant_classification': "5_prime_UTR"}, "5'UTR", inplace=True)
    mutation_merged.replace({'variant_classification': "Frameshift"}, "Frameshift_Variant", inplace=True)
    mutation_merged.replace({'variant_classification': "intergenic_variant"}, "Silent", inplace=True)

   # mutation_merged['variant_classification'] = mutation_merged['variant_classification'].replace("Nonsense", "Nonsense_Mutation", inplace=True)
    #mutation_merged['variant_classification'] = mutation_merged['variant_classification'].replace('intron', 'Intron', inplace=True)
    #mutation_merged['variant_classification'] = mutation_merged['variant_classification'].replace("synonymous", "Silent", inplace=True)
    #mutation_merged['variant_classification'] = mutation_merged['variant_classification'].replace("Inframe_Del", "In_Frame_Del", inplace=True)
    #mutation_merged['variant_classification'] = mutation_merged['variant_classification'].replace("5_prime_UTR", "5' UTR", inplace=True)
    #mutation_merged['variant_classification'] = mutation_merged['variant_classification'].replace("intergenic_variant", "Silent", inplace=True)
    mutation_merged_select = mutation_merged[['entrez_id', 'Sample_ID_Tumor', 'Name', 'variant_classification']]
    #merge with improve_ids 
    samples['other_id_no_dash'] = samples['other_id'].str.replace("-2", "_2")
    mutation_merged_2 = mutation_merged_select.merge(samples, left_on='Sample_ID_Tumor', right_on='other_id_no_dash', how='left')
    # select desired columns - entrez_id, improve_sample_id, mutation, variant_classificaton, source, study
    mutation_merged_2['other_id_source'] = "Synapse"
    mutation_merged_2['study'] = "Landscape of Sarcoma"
    mutationData = mutation_merged_2[['entrez_id',  'Name', 'variant_classification',  'improve_sample_id', 'study']]
    mutationData =mutationData.rename({"Name": "mutation"}, axis=1)
    # drop duplicates
    mutationData = mutationData.drop_duplicates()
    return mutationData


if __name__ == "__main__":
    print('in main')
    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of omics data files for the Sarcoma PDO project")
    parser.add_argument('-s', '--samples', help='Path to sample file',default=None)
    parser.add_argument('-g', '--genes', help='Path to genes file', default = None)
    parser.add_argument('-c', '--copy', help='Flag to capture copy number data', action='store_true', default=False)
    parser.add_argument('-m', '--mutation', help='Flag to capture mutation data', action='store_true', default=False)
    parser.add_argument('-e', '--expression', help='Flag to capture transcriptomic data', action='store_true', default=False)

    parser.add_argument('-t', '--token', help='Synapse token')

    args = parser.parse_args()
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)


    genes = pd.read_csv(args.genes)
    samples = pd.read_csv(args.samples)

    if args.expression:
        download_and_format_transcriptomic(synObject, genes, samples).to_csv("/tmp/sarcpdo_transcriptomics.csv", index=False)

   # if args.copy: 
    #    download_and_format_copy_number(synObject, genes, samples).to_csv('sarcpdo_copynumber.csv', index=False)

    if args.mutation:
        download_and_format_genomic_mutation(synObject, genes, samples).to_csv('/tmp/sarcpdo_mutation.csv', index=False)
    
          # validate with: linkml validate -s coderdata/schema/coderdata.yaml ~/Downloads/sarcpdo_samples.csv


    # command line testing: python3 01_createSarcPDOOmicsFiles.py -t $SYNAPSE_AUTH_TOKEN -s dev-environment/sarcpdo_samples.csv -g genes.csv -e
