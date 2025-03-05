import pandas as pd
import numpy as np
import wget
import os
import gzip
import requests
import argparse

def download_rnaseq(geo_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE65253&format=file&file=GSE65253%5Fcol%5Ftum%5Forg%5Fmerge%2Ecsv%2Egz", save_path = None): 
    """
    Retrieve data from a given GEO URL and identify the downloaded file by its name.

    This function uses the wget tool to download a file from the provided GEO URL.
    By comparing the directory contents before and after the download, 
    it identifies the newly downloaded file's name.

    Parameters
    ----------
    geo_url : str
        The GEO URL pointing to the data to be downloaded. Default is from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65253
    
    save_path : string
        Local path where the downloaded file will be saved.

    Returns
    -------
    None
    """
    
    response = requests.get(geo_url)
    with open(save_path, 'wb') as f:
        f.write(response.content)
    return

def download_sequencing_data(url = "https://www.cell.com/cms/10.1016/j.cell.2015.03.053/attachment/c2c327eb-2b79-493b-b658-9d45d8efdd74/mmc2.xlsx", save_path = None):
    """ 
    Download a file from a raw URL and save it to a local path.
    
    Parameters
    ----------
    raw_url : string
        The raw GitHub URL to download the file from.
        
    save_path : string
        Local path where the downloaded file will be saved.
        
    Returns
    -------
    None
    """
    
    response = requests.get(url)
    with open(save_path, 'wb') as f:
        f.write(response.content)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='###')
    
    parser.add_argument('-D', '--download',action='store_true', default=False, help='Download RNA seq and sequencing data from GEO and supplemental materials from https://www.cell.com/cell/fulltext/S0092-8674(15)00373-6#mmc2')
    parser.add_argument('-p', '--path', nargs='?',type=str, default='', const='', help='Local path to save downloaded RNA seq and sequencing data.')



    args = parser.parse_args()


    ###########################

    if args.download:
        if args.path is None or args.path=='':
            print("No path was provided to save data to. Cannot download data.")
            exit()
        else:
            print("Path to download data to was provided. Downloaded RNA seq and sequencing data to {}".format(args.path))
            # Download RNA seq data
            download_rnaseq(save_path = args.path)
            # Download sequencing data
            download_sequencing_data(save_path = args.path)
