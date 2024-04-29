'''
gets nci60 data from 10/2023 release

'''

import polars as pl
import argparse



conc_data='https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=11&modificationDate=1712351454136&api=v2'


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--sampleFile',dest='samplefile',default=None,help='DepMap sample file') 
    parser.add_argument('--drugfile',dest='dfile',default=None,help='Drug database')
    
    opts = parser.parse_args()
    
    samplefile = opts.samplefile
    drugfile = opts.dfile



if __name__=='__main__':
    main()
