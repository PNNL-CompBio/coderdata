
import os
import pandas as pd
import wget
import argparse
import synapseclient as sc
import math


def main():
    ##current AUC values are here: https://aacr.figshare.com/ndownloader/files/39996295 tabs 2 and 3
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pat', help='Synapse authentication token')
    parser.add_argument('-s', '--samples', help='Sample mapping file for panc pdo samples')
    parser.add_argument('-d', '--drugs', help='Drug mapping file for panc pdo samples')
    parser.add_argument('-o', '--output', default = '/tmp/pancpdo_doserep.tsv',help='Output file to be read into curve fitting code')

    args = parser.parse_args()
    newdata = get_data(args.pat)
    newdata = newdata.rename(columns={'Organoid':'other_id','Drug':'chem_name','Dose':'DOSE','FracResponse':'GROWTH','Passage':'time'})
    newdata = newdata[['other_id','chem_name','DOSE','GROWTH']]
    newdata[['time']]='120'
    newdata[['time_unit']]='hours'
    newdata[['study']]='pancpdo'
    newdata[['source']]='TiriacEtAl2018'
#    'source', 'improve_sample_id', 'Drug', 'study','time','time_unit'
    mappedresponse = map_drugs_to_samps(newddata,args.drugs,args,samples)
    mappedresponse.to_csv(args.output, sep='\t', index=False)

def map_to_drugs_samps(dose_rep,drugfile,sampfile):
    '''
    Collect dose response data frame, map drugs and organoids to improve drug and sample ids
    '''
    drugs = pd.read_csv(drugfile, sep='\t')
    samps = pd.read_csv(sampfile)

    merged = dose_rep.merge(drugs).merge(samps)
    merged = merged[['improve_sample_id','improve_drug_id','DOSE','GROWTH','time','time_unit','study','source']]
    merged = merged.rename(columns={'improve_drug_id':'Drug'})
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
    rtab = responses.melt(id_vars = responses.columns[0:4],value_vars=responses.columns[4:10], var_name='Drug',value_name='Response')
    
    ##rename the drugs
    rtab[['Drug','Rep']]=rtab['Drug'].str.lower().str.split('.',expand=True)
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
    rtab["MaxRep"] = rtab.groupby(['Drug','Organoid','Rep']).Response.transform('max')
    rtab['PercResponse'] = (rtab.Response/rtab.MaxRep)*100.00


    ##dosenum isa dummy value to use for merging since we need to repeat the concentrations over and over
    dosenum = [a for a in range(15)]
    rtab['Dosenum']=dosenum*int(rtab.shape[0]/15)
               
    ##merge the concentrations
    concs = concs.dropna().melt(value_vars=concs.columns,var_name='Drug',value_name='Dose')
    concs.Drug=concs.Drug.str.lower()
    concs['Dosenum'] = dosenum*int(concs.shape[0]/15)##creating dosenum here to merge

    
    return rtab.merge(concs)

if __name__=='__main__':
    main()
