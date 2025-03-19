import pandas as pd
import numpy as np

def parse_mmc2(mmc2_excel_path):
    """
    Creates sample file from sequencing data excel file. Checks the input sample file against previous sample files to make sure 
    there are no clashing sample names and assigns improved ID's starting from where previous sample sheet left off.
    
    Parameters
    ----------
    mmc2_excel_path : string
        Path to downloaded excel file. You can get this from the output of the function: download_sequencing_data()

    Returns
    -------
    mutations_data : pd.DataFrame
        A DataFrame containing mutations data.

    copy_number_data : pd.DataFrame
        A DataFrame containing copy number data.

    """
    mmc2_excel = pd.ExcelFile(open(mmc2_excel_path, 'rb'))
    mutations_data = pd.read_excel(mmc2_excel, 'TableS1J-Somatic mutations') # table with somatic mutation information 
    copy_number_data = pd.read_excel(mmc2_excel, 'TableS1D-Segmented_CN') # table with copy number information


    return(mutations_data, copy_number_data)

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

    
    # map columns in mutations data to their improved id
    mapped_mutation_data = mutations_data.loc[:,['Entrez_Gene_Id', 'Hugo_Symbol','Genome_Change','Variant_Classification','Tumor_Sample_Barcode']]
    mapped_mutation_data['Tumor_Sample_Barcode'] = mapped_mutation_data['Tumor_Sample_Barcode'].str.split('-',n = 1,expand=True).iloc[:,1]
    mapped_mutation_data = pd.merge(mapped_mutation_data, improve_id_data[['other_id','improve_sample_id']], how = 'left', left_on= "Tumor_Sample_Barcode", right_on= "other_id")
    
    # the data's variant classification matches scheme well, except "Non-coding_Transcript".  let's change those to RNA
    mapped_mutation_data.loc[mapped_mutation_data['Variant_Classification'] == "Non-coding_Transcript",'Variant_Classification'] = "RNA"

    # some of the entrezID's have ?, but it must be an integer.  check to see if we have an entrez ID using the hugo symbol and mapping using entrez_data. if there's any leftover, put in nothing
    questions = mapped_mutation_data[mapped_mutation_data['Entrez_Gene_Id'] == "?"].reset_index()
    fixed_entrez = pd.merge(questions, entrez_data[['entrez_id','other_id','gene_symbol']], how='inner', left_on="Hugo_Symbol", right_on="other_id")
    for index_val in fixed_entrez['index'].values:
        mapped_mutation_data.loc[index_val,'Entrez_Gene_Id'] = fixed_entrez[fixed_entrez['index'] == index_val]['entrez_id'].values
    mapped_mutation_data.loc[mapped_mutation_data['Entrez_Gene_Id'] == "?",'Entrez_Gene_Id'] = np.nan  # any leftover "?" change them to np.nan bc schema requires them to be int

    # clean up column names 
    mapped_mutation_data = mapped_mutation_data.rename(columns={'Entrez_Gene_Id':'entrez_id','Genome_Change':'mutation','Variant_Classification':'variant_classification'})
    mapped_mutation_data = mapped_mutation_data.drop(columns=['Hugo_Symbol','Tumor_Sample_Barcode','other_id'])
    mapped_mutation_data['source'] = "vandeWetering_2015"
    mapped_mutation_data['study'] = "CRC_Organoids"

    return(mapped_mutation_data)


def map_transcriptomics():
    return()


def map_copy_number():
    return()