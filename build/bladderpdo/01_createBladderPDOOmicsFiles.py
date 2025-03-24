import synapseclient
import pandas as pd
import numpy as np
import argparse
import os
import wget
import gzip
import subprocess
import math

def get_copy_call(a):
    """
    Heler Function - Determine copy call for a value.
    """

    if a is None:
        return float('nan')

    if math.isnan(a):
        return float('nan')

    a_val = math.log2(float(a)+0.000001)
    if a_val < 0.5210507:
        return 'deep del'
    elif a_val < 0.7311832:
        return 'het loss'
    elif a_val < 1.214125:
        return 'diploid'
    elif a_val < 1.422233:
        return 'gain'
    else:
        return 'amp'

    return pd.Series([get_copy_call(a) for a in arr])

def get_bladder_pdo_transcriptomics(GEO_id_link_table, samples, genes):

    bladderpdo_url ='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103990/suppl/GSE103990_Normalized_counts.txt.gz'
    transcriptomic_txt = wget.download(bladderpdo_url)
    transcriptomics = pd.read_csv(transcriptomic_txt, compression='gzip', sep="\t")
    subprocess.call (["/usr/local/bin/Rscript", "--vanilla", "obtainGSMidLink.R"])

    GEO_ids_link = pd.read_csv("./gsmlinkDf.csv")
    fpkm_totals = transcriptomics.iloc[:, 1:43].sum()
    transcriptomics.iloc[:, 1:43] = transcriptomics.iloc[:, 1:43].div(fpkm_totals).mul(1e6)
    transcriptomics['ensembl'] = transcriptomics['Unnamed: 0'].str.split("_", expand=True)[0]
    mapped_df = transcriptomics.merge(genes[['entrez_id', 'other_id']].drop_duplicates(), left_on='ensembl', right_on='other_id', how='left')
    # transform data to long format

    mapped_df.drop('other_id', axis=1)
    value_variables = transcriptomics.columns[transcriptomics.columns.str.contains("M")]
    melted_txomics = mapped_df.melt(id_vars = "entrez_id", value_vars = value_variables, var_name='sample_name')
    # use info from GEO to get Sample IDS
    txomics_with_GEOid = melted_txomics.merge(GEO_ids_link, how = 'left', left_on = "sample_name", right_on='RNAid')
    # use samplesheet to link sample_ids to improve ids
    txomics_with_GEOid['sampleid'] = txomics_with_GEOid['sampleid'].str.replace("org", "Organoid_")
    txomics_with_GEOid['sampleid'] = txomics_with_GEOid['sampleid'].str.replace("tumor", "Tumor")
    txomics_with_improveid = txomics_with_GEOid.merge(samples, left_on="sampleid", right_on="other_id", how="left")
    final_transcriptomics = txomics_with_improveid[['entrez_id', 'value', 'improve_sample_id']]
    final_transcriptomics['source'] = "Gene Expression Omnibus"
    final_transcriptomics['study'] = "Lee etal 2018 Bladder PDOs"
    final_transcriptomics.rename({'value' : 'transcriptomics' })
    # remove duplicates
    toreturn = final_transcriptomics.drop_duplicates()

    return toreturn

def get_bladder_pdo_mutations(synObject, samples, genes):
    print(samples.head)
    mutations = synObject.get("syn64765525")
    mutations_df = pd.read_csv(mutations.path, sep='\t')
    mutations_df['mutation'] = mutations_df['HGVSc'].str.split(":", expand=True)[1]
    #samplesheet = pd.read_csv(samples)
    selectioncols_mutations = mutations_df[['Entrez_Gene_Id',"Variant_Classification", "Tumor_Sample_Barcode", "mutation"]]
    merged_mutations = selectioncols_mutations.merge(samples, left_on="Tumor_Sample_Barcode", right_on="other_id", how="left")
    merged_mutations_renamed = merged_mutations.rename({"Entrez_Gene_Id" : 'entrez_id', 'Variant_Classification' : "variant_classification"}, axis=1)
    print(merged_mutations_renamed.head)
    final_mutations = merged_mutations_renamed[['entrez_id', "mutation", "variant_classification", "improve_sample_id"]]
    final_mutations['study'] = "Lee etal 2018 Bladder PDOs"
    print(final_mutations.head)
    return final_mutations

def get_bladder_pdo_copynumber(synObject, samples, genes):
    segfile = synObject.get("syn64765499")
    segfile_df = pd.read_csv(segfile.path, sep='\t')

    segfile_df.to_csv("bladder_segfile.csv")
    subprocess.call (["/usr/local/bin/Rscript", "--vanilla", "CNV-segfile-annotation.R", "bladder_segfile.csv", "bladder_annotated_segfile.csv"])
    copynumber = pd.read_csv("bladder_annotated_segfile.csv")
    copynumber['copy_number'] = np.exp2(copynumber['score'].div(2))*2
    copynumber['copy_call'] = [get_copy_call(a) for a in copynumber['copy_number']]
    copynumber_with_improveids = copynumber.merge(samples, left_on='ID', right_on = 'other_id', how='left')
    copynumber_with_correct_colnames = copynumber_with_improveids.rename({"ENTREZID":'entrez_id'}, axis=1)
    final_copynumber = copynumber_with_correct_colnames[['entrez_id', 'improve_sample_id', 'copy_number', 'copy_call']]
    final_copynumber['source'] = "Synapse"
    final_copynumber['study'] = "Lee etal 2018 Bladder PDOs"

    return final_copynumber




if __name__ == "__main__":
    print('in main')
    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of omics data files for the Bladder PDO project")
    parser.add_argument('-s', '--samples', help='Path to sample file',default=None)
    parser.add_argument('-g', '--genes', help='Path to genes file', default = None)
    parser.add_argument('-c', '--copy', help='Flag to capture copy number data', action='store_true', default=False)
    parser.add_argument('-m', '--mutation', help='Flag to capture mutation data', action='store_true', default=False)
    parser.add_argument('-e', '--expression', help='Flag to capture transcriptomic data', action='store_true', default=False)
    parser.add_argument('-i', '--geolink', help=".csv file that is the output of 'CNV-segfile-anotation.R")
    parser.add_argument('-t', '--token', help='Synapse token')

    args = parser.parse_args()
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)
    print('gene file is:')
    print(args.genes)
    print('sample file is :')
    print(args.samples)
    genes = pd.read_csv(args.genes)
    samples = pd.read_csv(args.samples)

    if args.expression:
        get_bladder_pdo_transcriptomics(args.geolink, samples, genes).to_csv("bladderpdo_transcriptomics.csv", index=False)

    if args.mutation:
        get_bladder_pdo_mutations(synObject, samples, genes).to_csv('bladderpdo_mutations.csv', index=False)
    
    if args.copy:
        get_bladder_pdo_copynumber(synObject, samples, genes).to_csv("bladderpdo_copynumber.csv", index=False)