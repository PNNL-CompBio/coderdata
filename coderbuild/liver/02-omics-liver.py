# load packages
import pandas as pd
import numpy as np
import os
import math
import argparse
import synapseclient
import mygene

# function for loading data
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
    proteomics_data = pd.read_excel(omics_excel, 'D. Proteomics')

    return(mutations_data, copynum_data, proteomics_data)


# function for getting rnaseq counts
def download_parse_rna_data(synID:str , save_path:str = None, synToken:str = None):
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
    rnaseq_data : pd.Dataframe
        A pandas dataframe containing rna sequencing counts data

    """
    
    syn = synapseclient.Synapse() 
    syn.login(authToken=synToken) 
    
    # Obtain a pointer and download the data 
    downloaded_data = syn.get(entity=synID, downloadLocation = save_path) 

    # Get the path to the local copy of the data file 
    rna_filepath = downloaded_data.path

    # Parse the downloaded excel file
    rna_excel = pd.ExcelFile(open(rna_filepath, 'rb'))
    rnaseq_data = pd.read_excel(rna_excel)

    return(rnaseq_data)



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
    mutation_data.columns = mutation_data.iloc[0]
    mutation_data = mutation_data.drop([0], axis=0)
    mutation_data['mutation']  = "g."+ mutation_data['Chromosome'] + ":" + np.where(mutation_data['Start_Position'] == mutation_data['End_Position'], mutation_data['Start_Position'].astype(str), mutation_data['Start_Position'].astype(str) + "_" + mutation_data['End_Position'].astype(str))
    for index, row in mutation_data.iterrows():
        if mutation_data.at[index,'Variant_Classification'].__contains__("Ins"):
            mutation_data.at[index,'mutation'] = mutation_data.at[index,'mutation'] + "ins" + mutation_data.at[index,'Tumor_Seq_Allele2']
        elif mutation_data.at[index,'Variant_Classification'].__contains__("Del"):
            mutation_data.at[index,'mutation'] = mutation_data.at[index,'mutation'] + "del" + mutation_data.at[index,'Tumor_Seq_Allele1']
        else:
            mutation_data.at[index,'mutation'] = mutation_data.at[index,'mutation'] + mutation_data.at[index,'Tumor_Seq_Allele1'] + ">" + mutation_data.at[index,'Tumor_Seq_Allele2']
            
    # map columns in mutations data to their improved id
    sample_mutation_data = pd.merge(mutation_data, improve_id_data[['other_id','improve_sample_id']], how='inner', left_on="Tumor_Sample_Barcode", right_on="other_id")

    # the data's variant classification matches scheme well, except "Non-coding_Transcript".  let's change those to RNA
    sample_entrez_mutation_data = pd.merge(sample_mutation_data, entrez_data[['entrez_id','other_id']], how='inner', left_on="Hugo_Symbol", right_on="other_id") # merge with our entrez database to see if we have additional matches

    # clean up column names and data types
    columns_to_drop = set(sample_entrez_mutation_data.columns) - set(['entrez_id','mutation','Variant_Classification','improve_sample_id'])
    mapped_mutation_data = sample_entrez_mutation_data.drop(columns=columns_to_drop)
    mapped_mutation_data = mapped_mutation_data.rename(columns={'Variant_Classification':'variant_classification'})
    mapped_mutation_data['source'] = "Synapse"
    mapped_mutation_data['study'] = "liver"
    mapped_mutation_data = mapped_mutation_data.dropna()
    mapped_mutation_data = mapped_mutation_data.astype({'entrez_id':'int','improve_sample_id':'int'})
    mapped_mutation_data = mapped_mutation_data.drop_duplicates()
    mapped_mutation_data = mapped_mutation_data[['entrez_id','mutation','variant_classification','improve_sample_id','study','source']]

    return(mapped_mutation_data)


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
    copy_number_data = copy_number_data.iloc[:,1:]
    copy_number_data.columns = copy_number_data.iloc[0]
    copy_number_data = copy_number_data.drop([0], axis=0)
    copynum_df = copy_number_data.drop(columns=['Hugo_Symbol','Cytoband'])


    # need to convert segment mean to copy number ratio, which is a 2 ^ x transformation, then melt df into 1 gene 1 sample per row
    long_cn_df = pd.melt(copynum_df, id_vars=['Gene ID'], value_vars=copynum_df.columns[copynum_df.columns != 'Gene ID'])

    # do copy_number calculation from score and get copy call column
    long_cn_df = long_cn_df.rename(columns={0:'other_id'})
    long_cn_df = long_cn_df.astype({'value':'float'})
    long_cn_df['copy_number'] = pow(2,long_cn_df['value'])*2
    long_cn_df['copy_call'] = [get_copy_call(a) for a in long_cn_df['copy_number']]

    # map ID to improve_ID
    improve_mapped_cn_df = pd.merge(long_cn_df, improve_id_data[['other_id','improve_sample_id']], how = 'left', on='other_id')

# clean up columns and data types
    improve_mapped_cn_df = improve_mapped_cn_df.drop(columns=['other_id','value'])
    improve_mapped_cn_df['source'] = "Synapse"
    improve_mapped_cn_df['study'] = "liver"
    improve_mapped_cn_df = improve_mapped_cn_df.dropna()
    improve_mapped_cn_df = improve_mapped_cn_df.rename(columns={'Gene ID':'entrez_id'})
    improve_mapped_cn_df = improve_mapped_cn_df.astype({'entrez_id':'int','improve_sample_id':'int'})
    improve_mapped_cn_df = improve_mapped_cn_df[['entrez_id','copy_number','copy_call','study','source','improve_sample_id']]
        
    return(improve_mapped_cn_df)

def map_transcriptomics(transciptomics_data, improve_id_data, entrez_data):

    # read in data
    if isinstance(transciptomics_data, pd.DataFrame) == False:
        transciptomics_data = pd.read_csv(transciptomics_data)

    if isinstance(improve_id_data, pd.DataFrame) == False:
        improve_id_data = pd.read_csv(improve_id_data)

    if isinstance(entrez_data, pd.DataFrame) == False:
        entrez_data = pd.read_csv(entrez_data)

    # first, convert genes, which are in ensembl id's to gene names
    transciptomics_data = transciptomics_data.rename(columns={'Unnamed: 0': 'stable_id'})
    mg = mygene.MyGeneInfo()
    ensembl_ids = transciptomics_data['stable_id'].values
    gene_info_list = mg.getgenes(ensembl_ids, fields='symbol')
    gene_df = pd.DataFrame.from_dict(gene_info_list)
    for_tpm = pd.merge(transciptomics_data, gene_df[['query','symbol']], how = 'inner', left_on= "stable_id", right_on= "query")
    for_tpm = for_tpm.dropna(subset=['symbol'])
    for_tpm = for_tpm.drop(columns=['query','stable_id'])
    for_tpm = for_tpm.rename(columns={'symbol':'stable_id'})
    for_tpm.to_csv("/tmp/counts_for_tpm_conversion.tsv", sep='\t')


    # run tpmFromCounts.py to convert counts to tpm
    os.system("python3 tpmFromCounts.py --counts /tmp/counts_for_tpm_conversion.tsv --genome_build https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.gtf.gz --gene_col stable_id --exclude_col stable_id --out_file /tmp/transcriptomics_tpm.tsv")

    # melt the df so there is one sample and gene per row
    long_transcriptomics_df = pd.read_csv("/tmp/transcriptomics_tpm.tsv",sep='\t')
    long_transcriptomics_df = pd.melt(long_transcriptomics_df, id_vars=['stable_id'], value_vars=long_transcriptomics_df.columns[long_transcriptomics_df.columns != 'stable_id'])
    long_transcriptomics_df = long_transcriptomics_df.rename(columns = {'value':'transcriptomics', 0:'sample_name'})


    # map gene names to entrez id's
    mapped_transcriptomics_df = pd.merge(long_transcriptomics_df, entrez_data[['other_id','entrez_id']].drop_duplicates(), how = 'inner', left_on= "stable_id", right_on= "other_id")
    mapped_transcriptomics_df = mapped_transcriptomics_df.dropna(subset=['entrez_id'])

    # mapping improve sample id'samples_df
    mapped_transcriptomics_df = pd.merge(mapped_transcriptomics_df, improve_id_data[['other_id','improve_sample_id']].drop_duplicates(), how = 'inner', left_on= "variable", right_on= "other_id")

    # clean up column names and data types
    mapped_transcriptomics_df = mapped_transcriptomics_df.drop(columns=['stable_id','variable','other_id_x','other_id_y'])
    mapped_transcriptomics_df['source'] = "Synapse"
    mapped_transcriptomics_df['study'] = "liver"
    mapped_transcriptomics_df = mapped_transcriptomics_df.dropna()
    mapped_transcriptomics_df = mapped_transcriptomics_df.astype({'entrez_id':'int','improve_sample_id':'int'})
    mapped_transcriptomics_df = mapped_transcriptomics_df[['entrez_id','transcriptomics','improve_sample_id','source','study']]

    return(mapped_transcriptomics_df)


def map_proteomics(proteomics_data, improve_id_data, entrez_data):

    # read in data
    if isinstance(proteomics_data, pd.DataFrame) == False:
        proteomics_data = pd.read_csv(proteomics_data, dtype=str, low_memory=False)

    if isinstance(improve_id_data, pd.DataFrame) == False:
        improve_id_data = pd.read_csv(improve_id_data)

    if isinstance(entrez_data, pd.DataFrame) == False:
        entrez_data = pd.read_csv(entrez_data)

    # first, replace colnames with first row and delete first row
    proteomics_data.columns = proteomics_data.iloc[0,:]
    proteomics_data = proteomics_data.iloc[1:]

    # melt the df so there is one sample and prot per row
    proteomics_data = proteomics_data.rename(columns = {proteomics_data.columns[0]:'gene_symbol'})
    long_prot_df = pd.melt(proteomics_data, id_vars=['gene_symbol'], value_vars=proteomics_data.columns[proteomics_data.columns != 'gene_symbol'])
    long_prot_df = long_prot_df.rename(columns = {0:'sample_name', 'value':'proteomics'})


    # Ensure both columns are string types for merging
    long_prot_df['gene_symbol'] = long_prot_df['gene_symbol'].astype(str)
    entrez_data['other_id'] = entrez_data['other_id'].astype(str)

    # map gene names to entrez id's
    mapped_proteomics_df = pd.merge(long_prot_df, entrez_data[['other_id','entrez_id']].drop_duplicates(), how = 'inner', left_on= "gene_symbol", right_on= "other_id")

    mapped_proteomics_df = mapped_proteomics_df.dropna(subset=['entrez_id'])

    # mapping improve sample id'samples_df
    mapped_proteomics_df = pd.merge(mapped_proteomics_df, improve_id_data[['other_id','improve_sample_id']].drop_duplicates(), how = 'inner', left_on= "sample_name", right_on= "other_id")

    # clean up column names and data types
    mapped_proteomics_df = mapped_proteomics_df.drop(columns=['gene_symbol','sample_name','other_id_x','other_id_y'])
    mapped_proteomics_df['source'] = "Synapse"
    mapped_proteomics_df['study'] = "liver"
    mapped_proteomics_df = mapped_proteomics_df.dropna()
    mapped_proteomics_df = mapped_proteomics_df.astype({'entrez_id':'int','improve_sample_id':'int'})
    mapped_proteomics_df = mapped_proteomics_df[['entrez_id','proteomics','improve_sample_id','source','study']]

    return(mapped_proteomics_df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='###')

    # arguments for file paths
    parser.add_argument('-g', '--genes', type=str, default=None, help='Path to genes.csv.  Can be obtained using this docker container: https://github.com/PNNL-CompBio/coderdata/blob/0225c52b861dcd6902521228731c54a61768bcd6/build/genes/README.md#L4')
    parser.add_argument('-i', '--ids', type=str, default=None, help='Path to sample Ids')
    parser.add_argument('-t', '--token', type=str, default=None, help='Synapse Token')


    # arguments for what data to process
    parser.add_argument('-P', '--parse', action = 'store_true', default=False, help='Parse excel file with data')
    parser.add_argument('-T', '--transcriptomics', action = 'store_true', default=False, help='Generate transcriptomics data')
    parser.add_argument('-M', '--mutations', action = 'store_true', default=False, help='Generate mutations data')
    parser.add_argument('-C', '--copy_number', action = 'store_true', default=False, help='Generate copy number data')
    parser.add_argument('-R', '--proteomics', action = 'store_true', default=False, help='Generate proteomics data')

    args = parser.parse_args()


    ###########################

    if args.parse:
        print("Parsing excel file.")
        # Download and parse rnaseq data
        rnaseq_df = download_parse_rna_data(synID="syn68327513", synToken = args.token, save_path="/tmp/")
        # Download rest of omics data
        mutation_data, copynum_df, proteomics_df= download_parse_omics_data(synID="syn66401303", synToken = args.token, save_path="/tmp/")
        # Save mutation and copy number data into csv format
        rnaseq_df.to_csv("/tmp/raw_rnaseq_data.csv")
        mutation_data.to_csv("/tmp/raw_mutation_data.csv")
        copynum_df.to_csv("/tmp/raw_copynum_data.csv")
        proteomics_df.to_csv("/tmp/raw_proteomics_data.csv")



    if args.transcriptomics:
        if args.genes is None or args.genes=='':
            print("No genes data provided. Exiting script.")
            exit()
        if args.ids is None or args.ids=='':
            print("No samples data provided. Exiting script.")
            exit()
        else:
            print("Starting transcriptomics data.")
            transcriptomics_df = map_transcriptomics(transciptomics_data = "/tmp/raw_rnaseq_data.csv", improve_id_data = "/tmp/liver_samples.csv", entrez_data = "/tmp/genes.csv")
            transcriptomics_df.drop_duplicates(inplace=True)
            transcriptomics_df.to_csv("/tmp/liver_transcriptomics.csv", index=False)
    
    if args.mutations:
        if args.genes is None or args.genes=='':
            print("No genes data provided. Exiting script.")
            exit()
        if args.ids is None or args.ids=='':
            print("No samples data provided. Exiting script.")
            exit()
        else:
            print("Starting mutations data.")
            mutation_df = map_mutations(mutation_data = "/tmp/raw_mutation_data.csv", improve_id_data = "/tmp/liver_samples.csv", entrez_data = "/tmp/genes.csv")
            mutation_df.to_csv("/tmp/liver_mutations.csv", index=False)
    
    if args.copy_number:
        if args.genes is None or args.genes=='':
            print("No genes data provided. Exiting script.")
            exit()
        if args.ids is None or args.ids=='':
            print("No samples data provided. Exiting script.")
            exit()
        else:
            print("Starting copy number data.")
            copy_number_df = map_copy_number(copy_number_data = "/tmp/raw_copynum_data.csv", improve_id_data = "/tmp/liver_samples.csv", entrez_data = "/tmp/genes.csv")
            copy_number_df.drop_duplicates(inplace=True)
            copy_number_df.to_csv("/tmp/liver_copy_number.csv", index=False)

    if args.proteomics:
        if args.genes is None or args.genes=='':
            print("No genes data provided. Exiting script.")
            exit()
        if args.ids is None or args.ids=='':
            print("No samples data provided. Exiting script.")
            exit()
        else:
            print("Starting proteomics data.")
            proteomics_df = map_proteomics(proteomics_data = "/tmp/raw_proteomics_data.csv", improve_id_data = "/tmp/liver_samples.csv", entrez_data = "/tmp/genes.csv")
            proteomics_df.to_csv("/tmp/liver_proteomics.csv", index=False)
    