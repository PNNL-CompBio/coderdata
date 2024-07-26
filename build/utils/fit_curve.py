#! /usr/bin/env python3

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager

import sys
import argparse
import numpy as np
import pandas as pd

from tqdm import tqdm
from itertools import islice
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
import multiprocessing

#import uno_data as ud

def format_coderd_schema(fname):
    """    formats output to comply with coderdata schema
    """
    df = pd.read_csv(fname,delimiter='\t')
    ##first rename Drug to improve_drug_id
    df2 = df.rename(columns={'Drug':'improve_drug_id'})
    new_df = pd.melt(df2,id_vars=['source','improve_sample_id','improve_drug_id','study','time','time_unit'],value_vars=['fit_auc','fit_ic50','fit_ec50','fit_r2','fit_ec50se','fit_einf','fit_hs','aac','auc','dss'],value_name='dose_response_value',var_name='dose_response_metric')

    new_df.to_csv(fname,sep='\t',index=False)

HS_BOUNDS_ORIG = ([0, 10**-12, 0], [1, 1, 4])

def hs_response_curve_original(x, einf, ec50, hs):
    """ from PharmacoDB supp. https://doi.org/10.1093/nar/gkx911
        bounds:
          einf: [0, 1]       # fraction of cells not susceptible to drug
          ec50: [10^-12, 1]  # concentration to have half target receptors bound: [1pM, 1M]
          hs:   [0, 4]       # hill slope binding cooperativity
    """
    return einf + (1 - einf) / (1 + np.power(x/ec50, hs))


HS_BOUNDS = ([0, 0, 0], [1, 12, 4]) 
HS_BOUNDS_NEG = ([0, -3,-1],[1,8,0]) ## made hill slope forced to be negative
def response_curve(x, einf, ec50, hs):
    """ transformed the original function with ec50 in -log10(M) instead of M
    """
    return einf + (1 - einf) / (1 + 10 ** ((ec50 - x) * hs))


def response_integral(x, einf, ec50, hs):
    return (1 - einf) * np.log10(1 + 10 ** ((ec50 - x) * hs)) / hs + x


def compute_area(x1, x2, einf, ec50, hs, mode='trapz'):
    popt = (einf, ec50, hs)
    if mode == 'trapz':
        # trapezoidal numerical integrationcollapse
        xx = np.linspace(x1, x2, 100)
        yy = response_curve(xx, *popt)
        area = np.trapz(yy, xx, dx=0.01)
    else:
        # the integral function can be expressed analytically
        # but sometimes less accurate due to float precision issues
        area = response_integral(x2, *popt) - response_integral(x1, *popt)
    return area



'''
added back this function as a spot check of data
'''
def fit_exp(df_exp, title=None, dmin=None, dmax=None, save=False):
    if save:
        font = {'family' : 'normal',
                # 'weight' : 'bold',
                'size'   : 14}
        matplotlib.rc('font', **font)
        plt.figure(figsize=(12, 6))

    print(df_exp)
    xdata = df_exp.DOSE.astype(float)
    ydata = df_exp.GROWTH.astype(float)
    # ydata = df_exp.GROWTH.clip(lower=0, upper=1.0).astype(float)

    # print(xdata)
    # print(ydata)

    popt, pcov = response_curve_fit(xdata, ydata)
    metrics = compute_fit_metrics(xdata, ydata, popt, pcov)

    if popt is None:
        return metrics

    dmin = dmin or xdata.min()
    dmax = dmax or xdata.max()
    xx = np.linspace(dmin, dmax, 100)
    yy = response_curve(xx, *popt)

    plt.xlim(dmax, dmin)
    plt.ylim(0, np.max([105, np.max(yy)]))
    plt.plot(xx, yy*100, 'r-', label='fit: einf=%.3f, ec50=%.3f, hs=%.3f' % tuple(popt))
    plt.plot(xdata, ydata.clip(lower=0, upper=1.0)*100, 'b*', label='')
    plt.xlabel('Dose (-log10(M))')
    plt.ylabel('Growth%')
    plt.title(title)
    plt.tight_layout()
    plt.legend()
    if save:
        plt.savefig('exp.png', dpi=360)
        plt.close()
    else:
        plt.show()

    return metrics.to_frame(name='metrics').T

def compute_fit_metrics(xdata, ydata, popt, pcov, d1=4, d2=10):
    '''
    xdata: dose data in log10(M)
    ydata: range from 0 to 1
    popt: fit curve metrics
    pcov: ??
    d1: minimum fixed dose in log10(M)
    d2: maximum fixed dose log10(M)
    '''
    if popt is None:
        cols = ['fit_auc','fit_ic50','fit_ec50','fit_ec50se','fit_r2','fit_einf','fit_hs','aac','auc','dss']#'auc ic50 ec50 ec50se R2fit rinf hs aac1 auc1 dss1'.split(' ')
        return pd.Series([np.nan] * len(cols), index=cols)
    einf, ec50, hs = popt
    perr = np.sqrt(np.diag(pcov))
    ec50se = perr[1]
    xmin = xdata.min()
    xmax = xdata.max()
    ypred = response_curve(xdata, *popt)
    r2 = r2_score(ydata, ypred)
    auc1 = compute_area(xmin, xmax, *popt) / (xmax - xmin)
    aac1 = 1 - auc1
    ic50 = ec50 - np.log10(0.5/(0.5-einf)) / hs if einf < 0.5 else np.nan
    ic90 = ec50 - np.log10(0.9/(0.1-einf)) / hs if einf < 0.1 else np.nan
    ic10 = ec50 - np.log10(0.1/(0.9-einf)) / hs if einf < 0.9 else np.nan
    ic10x = min(ic10, xmax)

    ##compute area under the ic10 to subtract from total
    int10x = compute_area(xmin, ic10x, *popt)
    ##old code - assumes a positive hill slope, otherwise doesn't seem to work.
    dss1 = (0.9 * (ic10x - xmin) - int10x) / (0.9 * (xmax - xmin)) if xmin < ic10x else 0
    auc = (response_integral(d2, *popt) - response_integral(d1, *popt)) / (d2 - d1)
    ##added by sara, i'm not sure where the above came from
    ## orig definition https://static-content.springer.com/esm/art%3A10.1038%2Fsrep05193/MediaObjects/41598_2014_BFsrep05193_MOESM1_ESM.pdf

    dss1 = (auc1-0.1*(ic10x-xmin)) / (0.9 * (xmax - xmin)) if xmin<ic10x else 0 #xmax > ic50 else 0
    dss2 = dss1/(1-einf) ##made this dss2 
    metrics = pd.Series({'fit_auc':auc, 'fit_ic50':ic50, 'fit_ec50':ec50,'fit_einf':einf,
                         'fit_ec50se':ec50se, 'fit_r2':r2, 'einf':einf, 'fit_hs':hs,
                         'aac':aac1, 'auc':auc1, 'dss':dss2}).round(4)
    return metrics
    ##and also this: https://github.com/bhklab/PharmacoGx/blob/master/R/computeDSS.R


def response_curve_fit(xdata, ydata, bounds=HS_BOUNDS_NEG):
    '''
     xdata: log10 molar concetnration
     ydata: value between 0 and 1 for response
     bounds: these are fixed in code, nto sure what they are for
    '''
    ydata = ydata.clip(lower=0, upper=1.0)
    popt, pcov = None, None
    nfev = 100 * 3
    while popt is None and nfev < 10000:
        # print(nfev)
        try:
            popt, pcov = curve_fit(response_curve, xdata, ydata, bounds=bounds, max_nfev=nfev)
            # popt, pcov = curve_fit(response_curve, xdata, ydata, bounds=bounds, max_nfev=nfev, method='dogbox')
        except RuntimeError:
            pass
        nfev *= 2
    return popt, pcov


def process_df(df, fname, sep='\t', ngroups=None):
    # df = df1.copy()
    i = 0
    header = None
    cols = ['source', 'improve_sample_id', 'Drug', 'study']
    groups = df.groupby(cols)
    f = open(fname, 'w')
    for name, group in tqdm(groups):
        # print(name)
        xdata = group.DOSE.astype(float)
        ##added the following 3 lines to acocunt for data normalized between 0 and 100 instead of 0 and 1
        ydata = group.GROWTH
      #  if max(ydata)>10:
      #      ydata = ydata/100.0
        ydata.clip(lower=0, upper=1.0).astype(float)
        popt, pcov = response_curve_fit(xdata, ydata)
        metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
        if header is None:
            header = cols + metrics.index.tolist()
            print(sep.join(header), file=f)
        print(sep.join(name), end=sep, file=f)
        print(sep.join([f'{x:.4g}' for x in metrics]), file=f)
        i += 1
        if ngroups and i >= ngroups:
            break
    f.close()


def process_single_drug(name_group_tuple):
    name, group = name_group_tuple
    xdata = group.DOSE.astype(float)
    ydata = group.GROWTH.clip(lower=0, upper=1.0).astype(float)
    popt, pcov = response_curve_fit(xdata, ydata)
    metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
    return name, metrics

def process_df_part(df, fname, beataml=False, sep='\t', start=0, count=None):
    cols = ['source', 'improve_sample_id', 'Drug', 'study','time','time_unit']
    groups = df.groupby(cols)
    count = count or (4484081 - start)
    groups = islice(groups, start, start+count)
    cores = multiprocessing.cpu_count()
    poolsize = round(cores-1)
    print('we have '+str(cores)+' cores and '+str(poolsize)+' threads')
    with multiprocessing.Pool(processes=poolsize) as pool:
        results = pool.map(process_single_drug, groups)

    with open(f'{fname}.{start}', 'w') as f:
        header = None
        for result in results:
            name, metrics = result
            if header is None:
                header = cols + metrics.index.tolist()
                print(sep.join(header), file=f)
            print(sep.join(str(n) for n in name), end=sep, file=f)
            print(sep.join(f'{x:.4g}' for x in metrics), file=f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='input file')
    parser.add_argument('--output', help='prefix of output file')
    parser.add_argument('--beataml', action='store_true', help='Include this if for BeatAML')

    args = parser.parse_args()
    print(args.input)
    df_all = pd.read_table(args.input)
    #drop nas
    df_all = df_all.dropna()
    ##pharmacoGX data is micromolar, we need log transformed molar
    df_all.DOSE = np.log10(df_all.DOSE*1000000)
    ##need data to be between 0 and 1, not 0 and 100
    df_all.GROWTH=df_all.GROWTH/100.00
    
    fname = args.output or 'combined_single_response_agg'
    process_df_part(df_all, fname, beataml=args.beataml)#, start=args.start, count=args.count)
    
#    if args.beataml == False:
    format_coderd_schema(fname+'.0')

if __name__ == '__main__':
    main()
