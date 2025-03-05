import pandas as pd
import numpy as np
import wget
import os
import gzip

def download_from_geo(url): 
    url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE65253&format=file&file=GSE65253%5Fcol%5Ftum%5Forg%5Fmerge%2Ecsv%2Egz"
    files_0 = os.listdir()
    wget.download(url)
    files_1 = os.listdir()
    new_file = str(next(iter(set(files_1) - set(files_0))))
    return new_file