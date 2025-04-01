import pandas as pd
import numpy as np
import os

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
    mutation_data : pd.DataFrame
        A DataFrame containing mutations data.

    copy_number_data : pd.DataFrame
        A DataFrame containing copy number data.

    """
    mmc2_excel = pd.ExcelFile(open(mmc2_excel_path, 'rb'))
    mutation_data = pd.read_excel(mmc2_excel, 'TableS1J-Somatic mutations') # table with somatic mutation information 
    copy_number_data = pd.read_excel(mmc2_excel, 'TableS1D-Segmented_CN') # table with copy number information


    return(mutation_data, copy_number_data)

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
    mapped_mutation_data = mutation_data.loc[:,['Entrez_Gene_Id', 'Hugo_Symbol','Genome_Change','Variant_Classification','Tumor_Sample_Barcode']]
    mapped_mutation_data['Tumor_Sample_Barcode'] = mapped_mutation_data['Tumor_Sample_Barcode'].str.split('-',n = 1,expand=True).iloc[:,1]
    mapped_mutation_data = pd.merge(mapped_mutation_data, improve_id_data[['other_id','improve_sample_id']], how = 'left', left_on= "Tumor_Sample_Barcode", right_on= "other_id")
    
    # the data's variant classification matches scheme well, except "Non-coding_Transcript".  let's change those to RNA
    mapped_mutation_data.loc[mapped_mutation_data['Variant_Classification'] == "Non-coding_Transcript",'Variant_Classification'] = "RNA"

    # some of the entrezID's have ?, but it must be an integer.  check to see if we have an entrez ID using the hugo symbol and mapping using entrez_data. if there's any leftover, put in nothing
    questions = mapped_mutation_data[mapped_mutation_data['Entrez_Gene_Id'] == "?"].reset_index()  # all rows with ? 
    fixed_entrez = pd.merge(questions, entrez_data[['entrez_id','other_id','gene_symbol']], how='inner', left_on="Hugo_Symbol", right_on="other_id") # merge with our entrez database to see if we have additional matches
    for index_val in fixed_entrez['index'].values:
        mapped_mutation_data.loc[index_val,'Entrez_Gene_Id'] = fixed_entrez[fixed_entrez['index'] == index_val]['entrez_id'].values  # for loop to replace these values with found entrez's
    mapped_mutation_data = mapped_mutation_data[mapped_mutation_data['Entrez_Gene_Id'] != "?"]  # remove rows with ? leftover

    # clean up column names and data types
    mapped_mutation_data = mapped_mutation_data.rename(columns={'Entrez_Gene_Id':'entrez_id','Genome_Change':'mutation','Variant_Classification':'variant_classification'})
    mapped_mutation_data = mapped_mutation_data.drop(columns=['Hugo_Symbol','Tumor_Sample_Barcode','other_id'])
    mapped_mutation_data['source'] = "vandeWetering_2015"
    mapped_mutation_data['study'] = "CRC_Organoids"
    mapped_mutation_data = mapped_mutation_data.astype({'entrez_id':'int'})

    return(mapped_mutation_data)


def map_transcriptomics(transciptomics_data, improve_id_data, entrez_data):

    # read in data
    if isinstance(transciptomics_data, pd.DataFrame) == False:
        transciptomics_data = pd.read_csv(transciptomics_data)

    if isinstance(improve_id_data, pd.DataFrame) == False:
        improve_id_data = pd.read_csv(improve_id_data)
    
    if isinstance(entrez_data, pd.DataFrame) == False:
        entrez_data = pd.read_csv(entrez_data)

    # move row names to a column called "stable_id" and format gene names to remove the chromosome num
    transciptomics_data['stable_id'] = transciptomics_data.index
    transciptomics_data['stable_id'] = transciptomics_data['stable_id'].str.split('__',n = 1,expand=True).iloc[:,0]
    transciptomics_data.to_csv("counts_for_tpm_conversion.csv")

    # run tpmFromCounts.py to convert counts to tpm
    os.system("python tpmFromCounts.py --counts counts_for_tpm_conversion.csv --genome_build https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.gtf.gz --gene_col stable_id --exclude_col stable_id --out_file transcriptomics_tpm.tsv")
    
    # get output from script (in tsv format) and average across organoids from each patient ]
    tpm_transciptomics_data = pd.read_csv("transcriptomics_tpm.tsv", sep="\t")
    tpm_transciptomics_data.index = tpm_transciptomics_data['stable_id']
    tpm_transciptomics_data = tpm_transciptomics_data.drop(columns=['stable_id'])
    transpose_transcriptomics = tpm_transciptomics_data.T
    transpose_transcriptomics['full_patient_id'] = transpose_transcriptomics.index
    transpose_transcriptomics['patient'] = transpose_transcriptomics['full_patient_id'].str.split('.',n = 1,expand=True).iloc[:,0]
    transpose_transcriptomics = transpose_transcriptomics.drop(columns=['full_patient_id'])
    transpose_transcriptomics = transpose_transcriptomics.groupby(by='patient').mean()

    # melt dataframe so that there is gene name and improve_sample_id per row
    transpose_transcriptomics = transpose_transcriptomics.reset_index()
    long_transcriptomics_df = pd.melt(transpose_transcriptomics, id_vars=['patient'], value_vars=transpose_transcriptomics.columns[transpose_transcriptomics.columns != 'patient'])
    long_transcriptomics_df = long_transcriptomics_df.rename(columns = {'value':'transcriptomics'})

    # map gene names to entrez id's 
    mapped_transcriptomics_df = pd.merge(long_transcriptomics_df, entrez_data[['other_id','entrez_id']].drop_duplicates(), how = 'left', left_on= "stable_id", right_on= "other_id")
    mapped_transcriptomics_df = mapped_transcriptomics_df.dropna(subset=['entrez_id'])

    # map patients to improve_sample_id 
    improve_id_tumor_org = improve_id_data[improve_id_data['other_id'].str.contains('Tumor-Organoid')] # mapping patient
    improve_id_tumor_org['patient'] = improve_id_tumor_org['other_id'].str.split('-',n = 1,expand=True).iloc[:,0].str.replace("T","").str.upper() # the way patient is written in transcriptomics data is slightly diff than that of our samples. edit to match. delete T's, uppercase A,B
    mapped_transcriptomics_df['patient'] = mapped_transcriptomics_df['patient'].str.upper() # make P in patients uppercase to merge
    mapped_transcriptomics_df = pd.merge(mapped_transcriptomics_df, improve_id_tumor_org[['patient','improve_sample_id']].drop_duplicates(), how = 'left', on='patient')
    
    # clean up column names and data types
    mapped_transcriptomics_df = mapped_transcriptomics_df.drop(columns=['stable_id','patient','other_id'])
    mapped_transcriptomics_df['source'] = "vandeWetering_2015"
    mapped_transcriptomics_df['study'] = "CRC_Organoids"
    mapped_transcriptomics_df = mapped_transcriptomics_df.astype({'entrez_id':'int','improve_sample_id':'int'})
    mapped_transcriptomics_df = mapped_transcriptomics_df[['entrez_id','transcriptomics','improve_sample_id','source','study']]

    return(mapped_transcriptomics_df)





def map_copy_number():
    return()