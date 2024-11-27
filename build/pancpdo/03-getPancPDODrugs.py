import pandas as pd
import os
import argparse




###figshare link:

filelink='https://aacr.figshare.com/ndownloader/files/39996295'
##get third tab and drugsa re listeda cross top

def main():
    parser = argparse.ArgumentParser(description='Download and match pancpdocdrugs')
    parser.add_argument('-d', '--prevDrugFile')
    parser.add_argument('-o', '--output', default = '/tmp/panpdc_drugs.tsv')

    

if __name__=='__main__':
    main()
