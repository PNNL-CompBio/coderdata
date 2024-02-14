#! /usr/bin/env python

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager
import pandas as pd

import sys
import argparse
import numpy as np
import pandas as pd

from tqdm import tqdm
from itertools import islice
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit

#import uno_data as ud

def format_coderd_schema(fname):
    '''
    formats output to comply with coderdata schema
    '''
    df = pd.read_table(fname)
    ##first rename Drug to improve_drug_id
    new_df = pd.melt(df,id_vars=['source','improve_sample_id','Drug','study','time','time_unit'],value_vars=['auc','ic50','ec50','ec50se','r2fit','hs','aac1','auc1','dss1'],value_name='dose_response_value',var_name='dose_response_metric')

    new_df.to_tsv(fname,sep='\t',index=False)

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
        cols = 'auc ic50 ec50 ec50se R2fit rinf hs aac1 auc1 dss1'.split(' ')
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
    int10x = compute_area(xmin, ic10x, *popt)
    dss1 = (0.9 * (ic10x - xmin) - int10x) / (0.9 * (xmax - xmin)) if xmin < ic10x else 0
    auc = (response_integral(d2, *popt) - response_integral(d1, *popt)) / (d2 - d1)
    metrics = pd.Series({'auc':auc, 'ic50':ic50, 'ec50':ec50,
                         'ec50se':ec50se, 'r2fit':r2, 'einf':einf, 'hs':hs,
                         'aac1':aac1, 'auc1':auc1, 'dss1':dss1}).round(4)
    return metrics


def response_curve_fit(xdata, ydata, bounds=HS_BOUNDS):
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


def fit_response(df_all, cell, drug, source, study=None, save=False):
#    cell_ids = ud.cell_name_to_ids(cell) or [cell]
#    drug_ids = ud.drug_name_to_ids(drug) or [drug]

    #df_exp = df_all[df_all.CELL.isin(cell_ids) & df_all.DRUG.isin(drug_ids)].copy()
    df_exp = df_all[(df_all.improve_sample_id == cell) & (df_all.Drug == drug)].copy()
    df_exp.GROWTH = (df_exp.GROWTH/2 + 0.5)
    df_exp = df_exp[df_exp.SOURCE == source]

    title = f'{cell} treated with {drug} in {source}'

    studies = df_exp.STUDY.unique()
    if len(studies) > 1:
        study = studies[study] if type(study) == int else study or studies[0]
        title += f' study {study}'
        df_exp = df_exp[df_exp.STUDY == study]

    return fit_exp(df_exp, title, save=save)


def show_dose_distribution(df_all):
    sources = df_all.SOURCE.unique()
    qs = [0, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.98, 1]
    series = []
    for src in sources:
        s = df_all[df_all.SOURCE == src].DOSE.quantile(qs)
        s.name = src
        series.append(s)
    df_dose = pd.concat(series, axis=1)
    return df_dose


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


def process_df_part(df, fname, sep='\t', start=0, count=None):
    header = None
    cols = ['source', 'improve_sample_id', 'Drug', 'study','time','time_unit']
    groups = df.groupby(cols)
    # count = count or (len(groups) - start)
    count = count or (4484081 - start)
    groups = islice(groups, start, start+count)
    f = open(f'{fname}.{start}', 'w')
    for name, group in tqdm(groups):
        #print(name)
        name = [str(n) for n in name]
        xdata = group.DOSE.astype(float)
#        ydata = group.GROWTH
#        ydata.clip(lower=0, upper=1.0).astype(float)
        ydata = group.GROWTH.clip(lower=0, upper=1.0).astype(float)
#        print(ydata)
        popt, pcov = response_curve_fit(xdata, ydata)
        metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
        if start == 0 and header is None:
            header = cols + metrics.index.tolist()
            print(sep.join(header), file=f)
        print(sep.join(name), end=sep, file=f)
        print(sep.join([f'{x:.4g}' for x in metrics]), file=f)
    f.close()





def process_chem_partner_data():
    df_cp = pd.read_csv('curve/ChemPartner_dose_response', sep='\t')
    df_cp = df_cp[df_cp.DRUG2.isnull() & df_cp.DOSE2.isnull()].drop(['DRUG2', 'DOSE2'], axis=1)
    df_cp = df_cp.rename(columns={'DRUG1':'DRUG', 'DOSE1':'DOSE'})
    df_cp.DOSE = -df_cp.DOSE
    # df_cp.GROWTH = df_cp.GROWTH/100
    df_cp.GROWTH = df_cp.GROWTH/200 + 0.5

    # process_df(df_cp, 'curve/ChemPartner_single_response_agg', ngroups=10)

    process_df(df_cp, 'curve/ChemPartner_single_response_agg.new')



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
    plt.plot(xx, yy*100, 'r-', label='fit: Einf=%.3f, EC50=%.3f, HS=%.3f' % tuple(popt))
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


def get_tableau20_colors():
    # tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
    #          (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
    #          (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
    #          (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
    #          (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                 (148, 103, 189), (44, 160, 44), (214, 39, 40), (255, 152, 150),
                 (152, 223, 138), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)
    return tableau20


def plot_curves(df_all, cell='LOXIMVI', drug='paclitaxel', study=None, max_reps=2, dmin=4, dmax=10, out=None):
#    cell_ids = ud.cell_name_to_ids(cell)
#    drug_ids = ud.drug_name_to_ids(drug)

    #df_exps = df_all[df_all.CELL.isin(cell_ids) & df_all.DRUG.isin(drug_ids)].copy()
    df_exps = df_all[(df_all['CELL']==cell) & (df_all['Drug']==drug)].copy()
    df_exps.GROWTH = (df_exps.GROWTH/2 + 0.5)

    title = f'{cell} treated with {drug}'
    out = out or f'{cell}-{drug}'

    # font = {'family': 'normal', 'size': 14}
    font = {'size': 14}
    matplotlib.rc('font', **font)
    plt.figure(figsize=(12, 6))
    colors = get_tableau20_colors()

    dmin = dmin or df_exps.DOSE.min()
    dmax = dmax or df_exps.DOSE.max()
    xx = np.linspace(dmin-0.1, dmax+0.1, 100)

    plt.xlim(dmax+0.1, dmin-0.1)
    plt.ylim(0, 105)
    plt.xlabel('Dose (-log10(M))')
    plt.ylabel('Growth%')
    plt.title(title)

    df_metrics = None
    rank = 0
    order = ['NCI60', 'CTRP', 'GDSC', 'CCLE', 'gCSI']
    sources = df_exps.SOURCE.unique().tolist() if study is None else study
    sources = sorted(sources, key=lambda x:order.index(x))

    for source in sources:
        studies = df_exps[df_exps.SOURCE == source].STUDY.unique()
        for i, study in enumerate(studies[:max_reps]):
            df_exp = df_exps[(df_exps.SOURCE == source) & (df_exps.STUDY == study)]
            xdata = df_exp.DOSE.astype(float)
            ydata = df_exp.GROWTH.astype(float)
            # ydata = df_exp.GROWTH.clip(lower=0, upper=1.0).astype(float)
            popt, pcov = response_curve_fit(xdata, ydata)
            metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
            if popt is None:
                continue
            color = colors[rank]
            rank = (rank + 1) % 20
            yy = response_curve(xx, *popt)
            label = source
            if len(studies) > 1:
                label += f' rep {i+1}'
            plt.plot(xx, yy*100, '-', color=color, label=label)
            plt.plot(xdata, ydata.clip(lower=0, upper=1.0)*100, '.', color=color, label='')
            if df_metrics is None:
                df_metrics = metrics.to_frame(name=label).T
            else:
                df_metrics = pd.concat([df_metrics, metrics.to_frame(name=label).T])

    plt.tight_layout()
    plt.legend()
    plt.savefig(f'{out}.png', dpi=360)
    plt.close()

    df_metrics.index.name = 'source'
    df_metrics.to_csv(f'{out}.csv', float_format='%.5g')
    print(f'Saved {out}.png and {out}.csv.')

    return df_metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='input file')
    parser.add_argument('--output', help='prefix of output file')

    args = parser.parse_args()
    print(args.input)
    df_all = pd.read_table(args.input)
    #drop nas
    df_all = df_all.dropna()
    #print(df_all)
    ##pharmacoGX data is micromolar, we need log transformed molar
    df_all.DOSE = np.log10(df_all.DOSE*1000000)
    ##need data to be between 0 and 1, not 0 and 100
    df_all.GROWTH=df_all.GROWTH/100.00
    
    fname = args.output or 'combined_single_response_agg'
    process_df_part(df_all, fname)#, start=args.start, count=args.count)
    format_coderd_schema(fname+'.0')

if __name__ == '__main__':
    main()
