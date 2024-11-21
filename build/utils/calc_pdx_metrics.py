'''
python script designed to compute pdx curve metrics from file
'''


from argparse import *
import pandas as pd
from tqdm import tqdm
import numpy as np

import statsmodels.api as sm
from statsmodels.formula.api import mixedlm


####mRECIST code
def tumor_volume_change(volume):
    """
    Calculate the percent change in tumor volume relative to the initial volume.
    
    Parameters:
    volume (np.ndarray or list): Array of tumor volumes.

    Returns:
    np.ndarray: Array of percent changes in tumor volume.
    """
    Vini = volume[0]
    return np.array([100 * (vt - Vini) / Vini for vt in volume])

def avg_response(volume_change):
    """
    Calculate the average response for each time point.
    
    Parameters:
    volume_change (np.ndarray or list): Array of volume percent changes.

    Returns:
    np.ndarray: Array of average responses.
    """
    cumsum = np.cumsum(volume_change)
    ar = cumsum / np.arange(1, len(volume_change) + 1)
    return ar

def check_numeric_int_char_zero(u):
    """
    Check if the input array is empty and return NA if true.

    Parameters:
    u (any): Input value or array.

    Returns:
    The input value or None if it is empty.
    """
    return None if len(u) == 0 else u

def get_best_response(time, response, min_time=None):
    """
    Determine the best response after a specified minimum time.

    Parameters:
    time (np.ndarray or list): Array of time points.
    response (np.ndarray or list): Array of response values.
    min_time (float or None): Minimum time after which tumor volume will be considered.

    Returns:
    dict: Dictionary containing the best response time, value, and index.
    """
    rtz = {"time": None, "value": None, "index": None}

    exdf = pd.DataFrame({"time": time, "response": response})
    exdf = exdf.dropna()
    exdf.reset_index()
    if min_time is not None:
        exdf_a = exdf[exdf['time'] >= min_time]

    if exdf.shape[0] == 0:
        return rtz
    exdf_a.reset_index()
 #   print(exdf.shape)
    min_indx = list(exdf.response).index(min(exdf_a.response))

    rtz['time'] = check_numeric_int_char_zero([list(exdf.time)[min_indx]])
    rtz['value'] = check_numeric_int_char_zero([list(exdf.response)[min_indx]])
    rtz['index'] = check_numeric_int_char_zero([min_indx])
    return rtz

def mrecist(time, volume, min_time=10, return_detail=False):
    """
    Compute the mRECIST for given volume response.

    Parameters:
    time (np.ndarray or list): Array of time points.
    volume (np.ndarray or list): Array of tumor volumes.
    min_time (float): Minimum time after which tumor volume will be considered.
    return_detail (bool): If True, return all intermediate values.

    Returns:
    str or dict: mRECIST classification or a dictionary with detailed intermediate values.
    """
    if volume[0] == 0:
        volume = np.array(volume) + 1

    exdf = {
        "volume_change": None,
        "average_response": None,
        "best_response": None,
        "best_response)time": None,
        "best_average_response": None,
        "best_average_response_time": None,
        "mRECIST": None
    }

    exdf["volume_change"] = tumor_volume_change(volume)
    exdf["average_response"] = avg_response(exdf["volume_change"])
    nxdf = {'metric':'mRESCIST','value':None}

    df = pd.DataFrame({"time": time, "volume": volume})
    df = df[df['time'] >= min_time]

    if df.shape[0] < 2:
        print(f"Warning: insufficient data after time {min_time}")
    else:
        br = get_best_response(time, exdf["volume_change"], min_time)
        exdf["best_response"] = br["value"]
        exdf["best_response_time"] = br["time"]

        bar = get_best_response(time, exdf["average_response"], min_time)
        exdf["best_average_response"] = bar["value"]
        exdf["best_average_response_time"] = bar["time"]

        best_response = exdf["best_response"][0]
        best_average_response = exdf["best_average_response"][0]

        mrecist = None
        nxdf = {'metric':'mRESCIST','value':mrecist}
#        exdf["mRECIST"] = mrecist

        if best_response is not None and best_average_response is not None:
            # Order of mRECIST assignment matters
            mrecist = "PD"

            if best_response < 35 and best_average_response < 30:
                mrecist = "SD"

            if best_response < -50 and best_average_response < -20:
                mrecist = "PR"

            if best_response < -95 and best_average_response < -40:
                mrecist = "CR"

#            exdf["mRECIST"] = mrecist
            nxdf = {'metric':'mRESCIST','value':mrecist}
    if not return_detail:
        return nxdf

    return nxdf

###AUC CODE

def trapz_auc(x, y):
    """
    Compute the area under the curve using trapezoidal rule.
    
    Parameters:
    x (np.ndarray): Array of x values.
    y (np.ndarray): Array of y values.
    
    Returns:
    float: Area under the curve.
    """
    n = np.arange(1, len(x))
    auc = np.dot((x[n] - x[n - 1]), (y[n] + y[n - 1])) / 2
    return auc

def AUC(time, volume, time_normalize=True):
    """
    Calculate the area under the curve (AUC).
    
    Parameters:
    time (np.ndarray): Array of time points recorded for the experiment.
    volume (np.ndarray): Array of volumes corresponding to the time points.
    time_normalize (bool): If True, AUC value will be divided by max time.
    
    Returns:
    dict: Dictionary containing the AUC value.
    """
    auc = trapz_auc(time, volume)
    #print(time)
    if time_normalize:
        auc = auc/np.max(time)
    return {"metric": "auc", "value": auc, 'time':np.max(time)}

def TGI(contr_volume, treat_volume,time):
    """
    Computes the tumor growth inhibition (TGI) between two time-volume curves.

    Parameters:
    contr_volume (np.ndarray or list): Volume vector for control.
    treat_volume (np.ndarray or list): Volume vector for treatment.

    Returns:
    dict: Dictionary containing the TGI value.
    """
    tgi = None
    if len(contr_volume) > 0 and len(treat_volume) > 0:
        tgi = (contr_volume[-1] - treat_volume[-1]) / (contr_volume[0] - treat_volume[0])
    
    # Simulated batch response class object
    rtx = {
        "metric": "TGI",
        "value": tgi,
        'time': np.max(time)
    }
    return rtx


def ABC(contr_time=None, contr_volume=None, treat_time=None, treat_volume=None):
    """
    Compute the area between two time-volume curves.
    
    Parameters:
    contr_time (np.ndarray): Array of time points for control.
    contr_volume (np.ndarray): Array of control volumes.
    treat_time (np.ndarray): Array of time points for treatment.
    treat_volume (np.ndarray): Array of treatment volumes.
    
    Returns:
    dict: Dictionary containing the area between curves.
    """
    con = {'name': 'auc', 'value': None}
    tre = {'name': 'auc', 'value': None}
    abc = None

    if contr_time is not None and contr_volume is not None:
        if len(contr_time) != len(contr_volume):
            raise ValueError("contr.time and contr.volume should have same length")
        con = AUC(contr_time, contr_volume)

    if treat_time is not None and treat_volume is not None:
        if len(treat_time) != len(treat_volume):
            raise ValueError("treat.time and treat.volume should have same length")
        tre = AUC(treat_time, treat_volume)

    abc = con['value'] - tre['value']
    return {"metric": "abc", "value": abc,'time':np.max(treat_time)}#, "control": con, "treatment": tre}


###LMM CODE
def lmm(time, volume, treatment, drug_name):
    """
    Compute the linear mixed model (lmm) statistics for a PDX batch.

    Parameters:
    data (pd.DataFrame): A DataFrame containing the batch data. Must contain the columns: model_id, volume, time, exp_type.

    Returns:
    dict: A dictionary containing the fit object and the specific coefficient value.
    """

    data = pd.DataFrame({'model_id':['model']*len(time),\
                         'volume':volume,\
                         'time':time,\
                          'exp_type':treatment})

    data = data.dropna()
                
    ##create data frame from these 4 vectors
    required_columns = ["model_id", "volume", "time", "exp_type"]
    
    if not all(column in data.columns for column in required_columns):
        raise ValueError("These columns must be present: 'model_id', 'volume', 'time', 'exp_type'")
    
    data['log_volume'] = np.log(data['volume'])
    
    # Define the formula for mixed linear model
    formula = 'log_volume ~ time*exp_type'
    
    # Fit the model
    model = mixedlm(formula, data, groups=data['model_id'])
    fit = model.fit()
    
    # Get the coefficient for the interaction term 'time:exp_type'
    #interaction_term = 'time:exp_type'
#    if interaction_term in fit.params:
#    time_coef_value = fit.params['time']
    #print(fit.params)
    i_coef_value = fit.params['time:exp_type[T.'+drug_name+']']
   # else:
   #     coef_value = None  # Handle the case when the interaction term is not present
    
    # Create a dictionary to store the fit object and the specific coefficient value
    result = {
        'fit': fit,
        'interaction_coef_value': i_coef_value
    }
    
    return {'metric': 'lmm','value': i_coef_value,'time':np.max(time)}

def main():
    parser=ArgumentParser()
    ###read in file with model id, volume, time, condition
    parser.add_argument('curvefile')
    parser.add_argument('--drugfile')
    parser.add_argument('--outprefix',default='/tmp/')
    
    args = parser.parse_args()
    
    tab = pd.read_csv(args.curvefile,sep='\t')
    drugs = pd.read_csv(args.drugfile,sep='\t')
    
    singles, combos = get_drug_stats(tab)

    ##join with drug ids
    expsing = singles.rename({'drug':'chem_name','metric':'drug_combination_metric','value':'drug_combination_value','sample':'improve_sample_id'},axis=1).merge(drugs,on='chem_name',how='left')[['improve_drug_id','improve_sample_id','drug_combination_metric','drug_combination_value']]
    expsing = expsing.dropna()

    combos[['drug1','drug2']]=combos.drug.str.split('+',expand=True)
    combos = combos.rename({'metric':'drug_combination_metric','value':'drug_combination_value','sample':'improve_sample_id'},axis=1).dropna()

    expcomb = combos.rename({'drug1':'chem_name'},axis=1).merge(drugs,on='chem_name',how='left').rename({'improve_drug_id':'improve_drug_1'},axis=1)[['improve_drug_1','drug2','improve_sample_id','drug_combination_metric','drug_combination_value']]
    expcomb = expcomb.rename({'drug2':'chem_name'},axis=1).merge(drugs,on='chem_name',how='left').rename({'improve_drug_id':'improve_drug_2'},axis=1)[['improve_drug_1','improve_drug_2','improve_sample_id','drug_combination_metric','drug_combination_value']]

    expcomb[['source']]='Synapse'
    expcomb[['study']]='MPNST PDX in vivo'

    expsing[['source']]='Synapse'
    expsing[['study']]='MPNST PDX in vivo'
    expsing.to_csv(args.outprefix+'_experiments.csv',index=False)
    expcomb.to_csv(args.outprefix+'_combinations.csv',index=False)
    

    
def get_drug_stats(df,control='control'):
    ##for each experiment, call group
    cols = ['experiment','model_id']
    groups = df.groupby(cols)
    singleres = []
    combores = []
    
    for name,group in tqdm(groups):
        #each group contains multiple treatments anda  control
        drugs = set(group.treatment)-set([control])
        print(name[0])
        print(drugs)
        mod = list(set(group.model_id))[0]
 #       print(set(group.model_id))
        ctl_data = group[group.treatment==control]
        ctl_time = np.array(ctl_data.time)
        ctl_volume = np.array(ctl_data.volume)

        ctl_auc = AUC(ctl_time,ctl_volume)
        for d in drugs:
            print(d)
            d_data = group[group.treatment==d]
            treat_time = np.array(d_data.time)
            treat_volume = np.array(d_data.volume)

            #get abc for group
            treat_auc = AUC(treat_time,treat_volume)
            treat_abc = ABC(ctl_time,ctl_volume,treat_time,treat_volume)
            #print(f"AUC: {treat_auc}")
            #print(f"ABC: {treat_abc}")
            treat_abc.update({'sample':mod,'drug':d,'time_unit':'days'})
            if '+' in d:
                combores.append(treat_abc)
            else:
                singleres.append(treat_abc)
            #lmm
            comb = pd.concat([ctl_data,d_data])
            lmm_res = lmm(comb.time, comb.volume, comb.treatment,d)
            lmm_res.update({'sample':mod,'drug':d,'time_unit':'days'})
            #print(f"LMM: {lmm_res}")
            if '+' in d:
                combores.append(lmm_res)
            else:
                singleres.append(lmm_res)

            #get tgi for group
            tg = TGI(ctl_volume,treat_volume,treat_time)
            tg.update({'sample':mod,'drug':d,'time_unit':'days'})
            #print(tg)
            if '+' in d:
                combores.append(tg)
            else:
                singleres.append(tg)

            
            #get mRECIST for group
            mr = mrecist(treat_time,treat_volume)
            mr.update({'sample':mod,'drug':d,'time_unit':'days'})
            if '+' in d:
                combores.append(mr)
            else:
                singleres.append(mr)

    sing = pd.DataFrame.from_records(singleres)
    comb = pd.DataFrame.from_records(combores)
    return sing,comb

if __name__=='__main__':
    main()
