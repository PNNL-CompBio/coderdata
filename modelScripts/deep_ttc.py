'''
deep_tty.py

Script to format data for the deepTTC algorithm
'''


import pandas as pd
import candle as cd


def get_data(type=['samples','drugs','gex','mut','experiment']):
    '''
    single function that stores the URL of the IMPROVE benchmark
    data, only used to test functionality
    '''
    url='https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/data/'
    if type=='samples':
        tab = pd.read_csv(url+'samples.csv')
    elif type=='drugs':
        tab = pd.read_csv(url+'drugs.tsv.gz',compression='gz',sep='\t')
    elif type=='gex':
        tab = pd.read_csv(url+'expression.csv.gz',comppression='gz',sep=',')
    elif type=='experiment':
        tab = pd.read_csv(url+'experiments.tsv.gz',compression='gz',sep=',')
    tab = pd.read_csv(url)


