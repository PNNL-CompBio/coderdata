'''
This script pulls down pre-computed curves and compares our fits with theirs
'''

import pandas as pd
import math
import argparse

# Figshare copy of the Tiriac et al. supplementary drug table (supp_4775187).
# The aacr.silverchair-cdn.com url previously here was a CloudFront SIGNED link
# whose Expires lapsed on 2025-01-27, returning 403 permanently; a fresh one
# expires about five weeks after issue. Figshare links do not expire.
# Tiriac H, et al. Cancer Discovery (2018) 8(9):1112-1129.
# doi:10.1158/2159-8290.CD-18-0349
tablink = 'https://ndownloader.figstatic.com/files/39996295'



def main():
    ##so far we have data for 'chemo' tab. how about the targeted tab?
    
    chemo = pd.read_excel(tablink,sheet_name=1)
    targeted = res = pd.read_excel(tablink,sheet_name=2)
    
    
    ##add in these scores to the drug file
    ##get drug file

    

if __name__=='__main__':
    main()
