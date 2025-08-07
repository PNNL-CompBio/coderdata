'''
This script pulls down pre-computed curves and compares our fits with theirs
'''

import pandas as pd
import math
import argparse

tablink = 'https://aacr.silverchair-cdn.com/aacr/content_public/journal/cancerdiscovery/8/9/10.1158_2159-8290.cd-18-0349/5/21598290cd180349-sup-199398_2_supp_4775187_p95dln.xlsx?Expires=1738004990&Signature=av8XadTm9AmI20O2Y7J7aHDtPbpluKJIfI5ubsoiYJ15D0zh5p1ltF4a7-DCSWTSMs-qX5TD09shxHeqkQ2NkLWHZsXoCD5KyREGhEgcDAvWZ1V9kwXDm0bjpINipAPPtC20oeuw6c~hPooF3Mtgzp4MzMCCjcVwfn05u27a0kS0yifBi11wQj3nmHlR3ym-2fYkFuqQtnNPCzH8-yIw21y0kTvXrNodAzC5pGA8qUK4PLxBt52xUIvTEPsPiPjXwBnDCfVsLGGdDYIY25lEPKiA403q6kFYvrSQ3bsTvM4kuvltb7yS4AXjK0-tthMOKbqq8~uREmJCcueADUF91g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'



def main():
    ##so far we have data for 'chemo' tab. how about the targeted tab?
    
    chemo = pd.read_excel(tablink,sheet_name=1)
    targeted = res = pd.read_excel(tablink,sheet_name=2)
    
    
    ##add in these scores to the drug file
    ##get drug file

    

if __name__=='__main__':
    main()
