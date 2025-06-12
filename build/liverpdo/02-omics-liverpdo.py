# load packages
import pandas as pd
import numpy as np
import os
import math
import argparse
import synapseclient

def download_parse_omics_data(synID:str , save_path:str = None, synToken:str = None):
    """ 
    Download omics data from Synapse at synapseID syn66401303. Requires a synapse token, which requires you to make a Synapse account
    and create a Personal Access Token.  More information here: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens 
    
    Parameters
    ----------
    synID : string
        SynapseID of dataset to download. Default is synapseID of the omics dataset.
        
    save_path : string
        atal path where the downloaded file will be saved.

    synToken : string
        Synapse Personal Access Token of user.  Requires a Synapse account. More information at: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
        
    Returns
    -------
    mutations_data : pd.Dataframe
        A pandas dataframe containing mutations data

    copynum_data : pd.Dataframe
        A pandas dataframe containing copy number data

    rnaseq_data : pd.Dataframe
        A pandas dataframe containing rna seq data

    proteomics_data : pd.Dataframe
        A pandas dataframe containing proteomics data
    """
    
    syn = synapseclient.Synapse() 
    syn.login(authToken=synToken) 
    
    # Obtain a pointer and download the data 
    downloaded_data = syn.get(entity=synID, downloadLocation = save_path) 

    # Get the path to the local copy of the data file 
    omics_filepath = downloaded_data.path

    # Parse the downloaded excel file
    omics_excel = pd.ExcelFile(open(omics_filepath, 'rb'))
    mutations_data = pd.read_excel(omics_excel, 'A. Mutation (LICOB)')
    copynum_data = pd.read_excel(omics_excel, 'B. CNV (LICOB)')
    rnaseq_data = pd.read_excel(omics_excel, 'C. RNA-seq (LICOB)')
    proteomics_data = pd.read_excel(omics_excel, 'D. Proteomics')

    return(mutations_data, copynum_data, rnaseq_data, proteomics_data)



def map_mutations(mutation_data, improve_id_data, entrez_data):
    """
    Maps mutation data to improved sample id's and entrez gene data. Also does some data formatting.
    
    Parameters
    ----------
    mutations_data : pd.Dataframe OR string
        Pandas dataframe object with mutation data OR path to csv with mutation data

    improve_id_data : pd.Dataframe OR string
        Pandas dataframe object with improve id data OR path to csv with improve id data.  This is one of the outputs of parse_mmc2()

    entrez_data : pd.Dataframe OR string
        Pandas dataframe object with entrez gene data OR path to csv with entrez gene data.  Use this code to get this file: https://github.com/PNNL-CompBio/coderdata/tree/e65634b99d060136190ec5fba0b7798f8d140dfb/build/genes 

    Returns
    -------
    mapped_mutation_data : pd.DataFrame
        A DataFrame containing the mapped mutations data with columns: entrez_id, mutation, variant_classification, improve_sample_id, source, study

    """

    # read in data
    if isinstance(mutation_data, pd.DataFrame) == False:
        mutation_data = pd.read_csv(mutation_data)

    if isinstance(improve_id_data, pd.DataFrame) == False:
        improve_id_data = pd.read_csv(improve_id_data)
    
    if isinstance(entrez_data, pd.DataFrame) == False:
        entrez_data = pd.read_csv(entrez_data)

    # create mutation names using chr, position, etc
    mutations_df.columns = mutations_df.iloc[0]
    mutations_df = mutations_df.drop([0], axis=0)
    mutations_df['mutation']  = "g."+ mutations_df['Chromosome'] + ":" + np.where(mutations_df['Start_Position'] == mutations_df['End_Position'], mutations_df['Start_Position'].astype(str), mutations_df['Start_Position'].astype(str) + "_" + mutations_df['End_Position'].astype(str))
    for index, row in mutations_df.iterrows():
        if mutations_df.at[index,'Variant_Classification'].__contains__("Ins"):
            mutations_df.at[index,'mutation'] = mutations_df.at[index,'mutation'] + "ins" + mutations_df.at[index,'Tumor_Seq_Allele2']
        elif mutations_df.at[index,'Variant_Classification'].__contains__("Del"):
            mutations_df.at[index,'mutation'] = mutations_df.at[index,'mutation'] + "del" + mutations_df.at[index,'Tumor_Seq_Allele1']
        else:
            mutations_df.at[index,'mutation'] = mutations_df.at[index,'mutation'] + mutations_df.at[index,'Tumor_Seq_Allele1'] + ">" + mutations_df.at[index,'Tumor_Seq_Allele2']
            
    # map columns in mutations data to their improved id
    sample_mutations_df = pd.merge(mutations_df, samples_df[['other_id','improve_sample_id']], how='inner', left_on="Tumor_Sample_Barcode", right_on="other_id")

    # the data's variant classification matches scheme well, except "Non-coding_Transcript".  let's change those to RNA
    sample_entrez_mutations_df = pd.merge(sample_mutations_df, entrez_df[['entrez_id','other_id']], how='left', left_on="Hugo_Symbol", right_on="other_id") # merge with our entrez database to see if we have additional matches

    # clean up column names and data types
    columns_to_drop = set(sample_entrez_mutations_df.columns) - set(['entrez_id','mutation','Variant_Classification','improve_sample_id'])
    mapped_mutations_df = sample_entrez_mutations_df.drop(columns=columns_to_drop)
    mapped_mutations_df = mapped_mutations_df.rename(columns={'Variant_Classification':'variant_classification'})
    mapped_mutations_df['source'] = "Synapse"
    mapped_mutations_df['study'] = "liverpdo"
    mapped_mutations_df = mapped_mutations_df.astype({'entrez_id':'int','improve_sample_id':'int'})
    mapped_mutations_df = mapped_mutations_df.drop_duplicates()
    mapped_mutations_df = mapped_mutations_df[['entrez_id','mutation','variant_classification','improve_sample_id','study','source']]

    return(mapped_mutations_df)


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


def map_copy_number(copy_number_data, improve_id_data, entrez_data):

    # read in data
    if isinstance(copy_number_data, pd.DataFrame) == False:
        copy_number_data = pd.read_csv(copy_number_data)

    if isinstance(improve_id_data, pd.DataFrame) == False:
        improve_id_data = pd.read_csv(improve_id_data)
    
    if isinstance(entrez_data, pd.DataFrame) == False:
        entrez_data = pd.read_csv(entrez_data)

    # get data ready 
    copy_number_data.columns = copy_number_data.iloc[0]
    copy_number_data = copy_number_data.drop([0], axis=0)
    copynum_df = copynum_df.drop(columns=['Hugo_Symbol','Cytoband'])


    # need to convert segment mean to copy number ratio, which is a 2 ^ x transformation, then melt df into 1 gene 1 sample per row
    long_transcriptomics_df = pd.melt(copynum_df, id_vars=['Gene ID'], value_vars=copynum_df.columns[copynum_df.columns != 'Gene ID'])

    # do copy_number calculation from score and get copy call column
    long_transcriptomics_df = long_transcriptomics_df.rename(columns={0:'other_id'})
    long_transcriptomics_df['copy_number'] = pow(2,long_transcriptomics_df['value'])*2
    long_transcriptomics_df['copy_call'] = [get_copy_call(a) for a in long_transcriptomics_df['copy_number']]

    # map ID to improve_ID
    improve_mapped_cn_df = pd.merge(long_transcriptomics_df, improve_id_data[['other_id','improve_sample_id']], how = 'left', on='other_id')

# clean up columns and data types
    improve_mapped_cn_df = improve_mapped_cn_df.drop(columns=['other_id','value'])
    improve_mapped_cn_df['source'] = "Synapse"
    improve_mapped_cn_df['study'] = "liverpdo"
    improve_mapped_cn_df = improve_mapped_cn_df.rename(columns={'Gene ID':'entrez_id'})
    improve_mapped_cn_df = improve_mapped_cn_df.astype({'entrez_id':'int','improve_sample_id':'int'})
    improve_mapped_cn_df = improve_mapped_cn_df[['entrez_id','copy_number','copy_call','study','source','improve_sample_id']]
        
    return(improve_mapped_cn_df)
