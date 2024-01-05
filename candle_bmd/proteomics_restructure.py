import pandas as pd


def import_github_data(url):
    """
    Import Data from github.

    This function simply wraps the pandas `read_csv` function.

    Parameters
    ----------
    url : string
        URL to pull data from.
    
    Returns
    -------
    pandas dataframe
        
    Notes
    --------
    If this fails, change VPN or network.
    """
    return pd.read_csv(url)
    
def reformat_gene_data(gene_data):
    """
    Reformat Gene Data

    This function modifies the index and changes the 'other_id' column to 'ENSG_col'.

    Parameters
    ----------
    gene_data : pandas dataframe
        dataframe to modify
    
    Returns
    -------
    updated pandas dataframe
     """
        
    gene_data.rename(columns={"other_id":"ENSG_col"},inplace = True)
    gene_data = gene_data.groupby('ENSG_col').first().reset_index()
    gene_data = gene_data.dropna(subset=['ENSG_col'])
    return gene_data
    
def merge_all_data(proteomics_data, sample_data, gene_data):
    """
    Merge all data into a single dataframe

    This function merges proteomics, sample, and gene data. It aggregates duplicate samples/proteins by mean. 
    Last it drops duplicates and NAs present in the other_id column.

    Parameters
    ----------
    proteomics_data : pandas dataframe
        
    sample_data : pandas dataframe
        
    gene_data : pandas dataframe
    
    Returns
    -------
    A single merged pandas dataframe
    """
        
    merged_data = proteomics_data.merge(sample_data[['improve_sample_id', 'other_id']], on='improve_sample_id', how='left')
    merged_data = merged_data.merge(gene_data[["entrez_id","gene_symbol","ENSG_col"]], on= "entrez_id",how = "left")

    # aggregate data by mean to join together duplicates
    merged_data = merged_data.groupby(['improve_sample_id', 'entrez_id',"gene_symbol","ENSG_col"]).agg({
        'proteomics': 'mean',
        'other_id': 'first'
    }).reset_index()
    merged_data.drop(columns=["improve_sample_id"],inplace=True)
    merged_data.drop_duplicates(inplace=True)
    merged_data = merged_data.dropna(subset=['other_id'])
    return merged_data

def restructure_data(merged_data):
    """
    Restructures merged dataframe to match desired output.

    Desired format has a 3 part header using ENSG identifier, gene symbol, and entrez id.
    
    Parameters
    ----------
    merged_data : pandas dataframe
    
    Returns
    -------
    A single restructured pandas dataframe
    """
        
    # Pivot merged dataframe
    restructured_data = merged_data.pivot(index='other_id', columns=["ENSG_col",'entrez_id',"gene_symbol"], values='proteomics').reset_index()

    # Adjust headers
    cols = list(restructured_data.columns)
    cols[0] = ('', '', '')
    header = pd.MultiIndex.from_tuples(cols)
    restructured_data.columns = header

    # Remove Duplicates by header
    restructured_data = restructured_data.T
    restructured_data = restructured_data.reset_index().drop_duplicates(subset=['level_1', 'level_2']).set_index(['level_0', 'level_1', 'level_2'])
    restructured_data = restructured_data.T

    # Adjust headers again
    cols = list(restructured_data.columns)
    cols[0] = ('', '', '')
    header = pd.MultiIndex.from_tuples(cols)
    restructured_data.columns = header
    return restructured_data


def write_results(data):
    """
    Write two TSVs as output. One output includes NA values. The second output fills NA values with zeros.

    This function wraps the pandas "to_csv" function.
    
    Parameters
    ----------
    data : pandas dataframe
    
    Returns
    -------
    None
    """
        
    data.to_csv("proteomics_restructure.tsv", index=False, sep="\t")
    data = data.fillna(0) 
    data.to_csv("proteomics_restructure_NA_as_zero.tsv", index=False, sep="\t")
    return

if __name__ == "__main__":
    
    # define urls
    proteomics_url = 'https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/cell_line/proteomics.csv.gz'
    samples_url = 'https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/cell_line/samples.csv'
    gene_url = 'https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/cell_line/genes.csv'

    # import data from github
    p_data =import_github_data(proteomics_url)
    s_data = import_github_data(samples_url)
    g_data = import_github_data(gene_url)

    # reformat gene data
    g_data = reformat_gene_data(g_data)

    # merge data
    merged_data = merge_all_data(p_data, s_data, g_data)

    # restructure data to match desired format
    restructured_data = restructure_data(merged_data)
    
    # write to csv 
    write_results(restructured_data)
