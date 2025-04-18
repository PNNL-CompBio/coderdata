
import os
import pandas as pd
import wget
import argparse
import synapseclient as sc
import math
import re

def main():
    ##current AUC values are here: https://aacr.figshare.com/ndownloader/files/39996295 tabs 2 and 3
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pat', help='Synapse authentication token')
    parser.add_argument('-s', '--samples', help='Sample mapping file for panc pdo samples')
    parser.add_argument('-d', '--drugs', help='Drug mapping file for panc pdo samples')
    parser.add_argument('-o', '--output', default = '/tmp/pancpdo_doserep.tsv',help='Output file to be read into curve fitting code')

    args = parser.parse_args()
    newdata = get_data(args.pat)
    newdata = newdata.rename(columns={'Organoid':'other_id','Drug':'chem_name','Dose':'DOSE','PercResponse':'GROWTH','Passage':'time'})
#    print(newdata)
    newdata = newdata[['other_id','chem_name','DOSE','GROWTH']]
    newdata[['time']]='120'
    newdata[['time_unit']]='hours'
    newdata[['study']]='pancpdo'
    newdata[['source']]='TiriacEtAl2018'
    print('collected doses and response for '+str(len(set(newdata.chem_name)))+' drugs and '+str(len(set(newdata.other_id)))+' samples')
#    'source', 'improve_sample_id', 'Drug', 'study','time','time_unit'
    mappedresponse = map_to_drugs_samps(newdata,args.drugs,args.samples)
    print('mapped doses and response for '+str(len(set(mappedresponse.Drug)))+' drugs and '+str(len(set(mappedresponse.improve_sample_id)))+' samples')
    mappedresponse.to_csv(args.output, sep='\t', index=False)

def map_to_drugs_samps(dose_rep,drugfile,sampfile):
    '''
    Collect dose response data frame, map drugs and organoids to improve drug and sample ids
    '''
    drugs = pd.read_csv(drugfile, sep='\t')
    samps = pd.read_csv(sampfile)

    merged = dose_rep.merge(drugs).merge(samps)

    merged = merged.rename(columns={'improve_drug_id':'Drug'}) 
    merged = merged[['improve_sample_id','Drug','DOSE','GROWTH','time','time_unit','study','source']].drop_duplicates()
    print(merged)
    return merged


def get_data(token):
    synid = 'syn64333325'

    syn = sc.login(authToken=token)
    fpath = syn.get(synid).path
    print(fpath)
    concs = pd.read_excel(fpath,sheet_name='concentrations')

    responses = pd.read_excel(fpath,sheet_name='Sheet1').dropna(axis=0,how='all')

    ##kludgy way of fixing rows so that all data is in each row
    newrows=[]
    org=''
    passage=''
    date=''
    pate=''
    responses = responses.fillna('').reset_index(drop=True)
    for rownum, row in responses.iterrows():
        if row['Organoid']!="":
            org = row['Organoid']
            passage = row['Passage']
            date = row['Date']
            pate = row['pate']
        newrows.append({'Organoid':org,'Passage':passage,'Date':date,'pate':pate})
        
    releft = pd.DataFrame(newrows)
    responses.Organoid = releft.Organoid
    responses.Passage = releft.Passage
    responses.pate = releft.pate
    responses.Date = releft.Date

    
    ##now melt the data into single columns
    rtab = responses.melt(id_vars = responses.columns[0:4],value_vars=responses.columns[4:20], var_name='Drug',value_name='Response')
    print('Collected results from '+str(len(set(rtab.Drug)))+' drugs and '+str(len(set(rtab.Organoid)))+' organoids')
    #print(set(rtab.Drug))
    ##rename the drugs
    rtab[['Drug','Rep']]=rtab['Drug'].str.lower().str.split('.',expand=True)
    rtab.Drug=[re.sub('-','',a) for a in rtab.Drug]
    #print(set(rtab.Drug))
    newrep=[]
    for r in rtab.Rep:
        if r is None:
            newrep.append(0)
        else:
            newrep.append(r)
    rtab.Rep=newrep

    ##renormalize values to max
    ##IMPORTANT: this is how we normalize without DMSO. We need to consider how we're doing this for EACH ORGANOID
    ##currently we take the max value of each orgnaoid/replicate.
    ##UPDATE: see belo
#    rtab["MaxRep"] = rtab.groupby(['Drug','Organoid','Rep']).Response.transform('max')
#    rtab['PercResponse'] = (rtab.Response/rtab.MaxRep)*100.00


    ##dosenum isa dummy value to use for merging since we need to repeat the concentrations over and over
    dosenum = [a for a in range(15)]
    rtab['Dosenum']=dosenum*int(rtab.shape[0]/15)

    ##The last dose (dosenum ==14) is the control value per Herve. we now must normalize to that

    dmso_vals = rtab[rtab.Dosenum==14][['Organoid','Drug','Rep','Response']].rename({'Response':'DMSO'},axis=1)
    full_res = rtab.merge(dmso_vals,on=['Organoid','Drug','Rep'])
    full_res['PercResponse'] = 100*(full_res.Response/full_res.DMSO)
    
    #print(set(rtab.Drug))
    ##merge the concentrations
    concs = concs.dropna().melt(value_vars=concs.columns,var_name='Drug',value_name='Dose')
    print(concs)
    concs.Dose = [d*10.0**6.0 for d in concs.Dose] ## convert M to uM here
    
    concs.Drug=concs.Drug.str.lower()
    concs['Dosenum'] = dosenum*int(concs.shape[0]/15)##creating dosenum here to merge
    #print(set(concs.Drug))
    
    return full_res.merge(concs)

if __name__=='__main__':
    main()
