import pandas as pd
import os
import argparse
import synapseclient



###figshare link:

filelink='https://aacr.figshare.com/ndownloader/files/39996295'
synid = 'syn64333325'
##get third tab and drugsa re listeda cross top


    

def retrieve_figshare_data(url):
    """
    Download data from a given Figshare URL.
    
    Parameters
    ----------
    url : string
        The Figshare URL to download data from.
    
    Returns
    -------
    string
        Name of the downloaded file.
    """
    
    files_0 = os.listdir()
    wget.download(url)
    files_1 = os.listdir()
    new_file = str(next(iter(set(files_1) - set(files_0))))
    return new_file

def main():
    parser = argparse.ArgumentParser(description='Download and match pancpdocdrugs')
    parser.add_argument('-d', '--prevDrugFile')
    parser.add_argument('-o', '--output', default = '/tmp/panpdc_drugs.tsv')

    auc_file = retrieve_figshare_data(filelink)

    tab = pd.read_excel(auc_file,sheet='')
    
    
if __name__=='__main__':
    main()
