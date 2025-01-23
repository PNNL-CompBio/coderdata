import pandas as pd
import synapseclient
import argparse
import math


def get_copy_call(a):
    """
    Helper Function - Determine copy call for a value.
    """

    if a is None:
        return float('nan')
    
    if math.isnan(a):
        return float('nan')
    
    a_val = a##math.log2(float(a)+0.000001) ###this should not be exponent, should be log!!! 2**float(a)
    if a_val < 0.0: #0.5210507:
        return 'deep del'
    elif a_val < 0.7311832:
        return 'het loss'
    elif a_val < 1.214125:
        return 'diploid'
    elif a_val < 1.731183:
        return 'gain'
    else:
        return 'amp'
    
    return pl.Series([get_copy_call(a) for a in arr])

def parseCNVFile(fpath, sampid, genes):
    log2data = pd.read_csv(fpath, sep='\t', header=None)
    log2data.columns = ['gene_symbol','copy_number','Region','Type','Pos']
    log2data['improve_sample_id']=sampid
    newdat =  pd.merge(log2data,genes)[['improve_sample_id','entrez_id','copy_number']].drop_duplicates()
    newdat['study']='pancpdo'
    newdat['source']='TiriacEtal'
    newdat = newdat[['improve_sample_id','entrez_id','copy_number','source','study']]
    newdat['copy_call'] = [get_copy_call(a) for a in newdat['copy_number']]
    return newdat


mutmap = {'CODON_CHANGE_PLUS_CODON_DELETION':'In_Frame_Del', ##this isn't a great mapping
          'CODON_CHANGE_PLUS_CODON_INSERTION':'In_Frame_Ins', ##this isn't a great mapping
          'CODON_DELETION':'In_Frame_Del',
          'CODON_INSERTION':'In_Frame_Ins',
          'DOWNSTREAM':"3'Flank",
          'FRAME_SHIFT':'Frameshift_Variant',
          'FRAME_SHIFT+SPLICE_SITE_ACCEPTOR+SPLICE_SITE_REGION+INTRON':'Frameshift_Variant',
          'FRAME_SHIFT+SPLICE_SITE_REGION':'Frameshift_Variant',
          'INTERGENIC':'IGR',
          'INTRON':'Intron',
          'NON_SYNONYMOUS_CODING':'Missense_Mutation',
          'NON_SYNONYMOUS_CODING+SPLICE_SITE_REGION':'Missense_Mutation',
          'SPLICE_SITE_ACCEPTOR+INTRON':'Splice_Site',
          'SPLICE_SITE_DONOR+INTRON':'Splice_Site',
          'SPLICE_SITE_REGION+INTRON':'Splice_Site',
          'SPLICE_SITE_REGION+NON_CODING_EXON_VARIANT':'Splice_Site',
          'SPLICE_SITE_REGION+SYNONYMOUS_CODING':'Silent',
          'START_GAINED+UTR_5_PRIME':'Start_Codon_Ins',
          'STOP_GAINED':'Stop_Codon_Ins',
          'STOP_GAINED+CODON_CHANGE_PLUS_CODON_INSERTION':'Stop_Codon_Ins',
          'SYNONYMOUS_CODING':'Silent',
          'UPSTREAM':"5'Flank",
          'UTR_3_PRIME':"3'UTR",
          'UTR_5_PRIME':"5'UTR"
          }

def parseMutFile(fpath, sampid,genes):
    '''
    move mutations to following headers:
     entrez_id, improve_sample_id, source, study, mutation, variant_classification
    '''
    mutfile = pd.read_csv(fpath,sep='\t')[['SNPEFF_GENE_NAME','SNPEFF_EFFECT','SNPEFF_CDS_CHANGE']]
    mutfile = mutfile.dropna(subset='SNPEFF_CDS_CHANGE')
    mutfile.columns  = ['gene_symbol','SNPEFF_EFFECT','mutation']
    fullfile = pd.merge(mutfile,pd.DataFrame({'SNPEFF_EFFECT':mutmap.keys(),'variant_classification':mutmap.values()}))
    fullfile = pd.merge(fullfile,genes)
    fullfile['improve_sample_id'] = sampid
    fullfile['source']='TiriacEtAl'
    fullfile['study']='pancpdo'
    fullfile = fullfile[['improve_sample_id','entrez_id','source','study','mutation','variant_classification']]
    fullfile = fullfile.dropna().drop_duplicates()
    return fullfile

def main():
    parser = argparse.ArgumentParser(description = 'Script that collects WES and CNV data from Synapse for Coderdata')
    parser.add_argument('-s', '--samples', help='Path to sample file',default=None)
    parser.add_argument('-g', '--genes', help='Path to genes file', default = None)
    parser.add_argument('-c', '--copy', help='Flag to capture copy number data', action='store_true', default=False)
    parser.add_argument('-m', '--mutation', help='Flag to capture mutation data', action='store_true', default=False)
    parser.add_argument('-t', '--token', help='Synapse token')

    args = parser.parse_args()
    if args.samples is None or args.genes is None:
        print('We need at least a genes and samples file to continue')
        exit()
    samps = pd.read_csv(args.samples)
    genes = pd.read_csv(args.genes)

    sc = synapseclient.login(args.token)
    ##to double check identifiers, we use transcriptomics data since that determines what samples were sequenced
    trans = pd.read_csv('/tmp/pancpdo_transcriptomics.csv.gz')
    tsamps = samps[samps.improve_sample_id.isin(trans.improve_sample_id)]
    print(samps.shape)
    print(tsamps.shape)
    
          
    missingsamples = []
    if args.copy:
        ##query synapse view for files
        cnvs = sc.tableQuery("select * from syn64608378 where parentId='syn64608163'").asDataFrame()
        alldats = []
        ##go through table and get every file
        for index,row in cnvs.iterrows():
            sid = row.id
            sname = row['name'].split('--')[0]
            print(sid,sname)
            path = sc.get(sid).path
            if sname in set(tsamps.other_id):
                print(sname+' in transcriptomics, using that id')
                sampid = tsamps.loc[tsamps.other_id==sname]['improve_sample_id'].values[0]
                missingsamples.append('copy,trans,'+sname)
            elif sname in set(samps.other_id):
                print(sname+' in samples but not transcriptomics, using other id')
                sampid = samps.loc[samps.other_id==sname]['improve_sample_id'].values[0]
                missingsamples.append("copy,notrans,"+sname)
            else:
                print('Missing sample id for '+sname,' skipping for now')
                missingsamples.append('copy,missed,'+sname)
                continue
            sampid = samps.loc[samps.other_id==sname]['improve_sample_id'].values[0]
            res = parseCNVFile(path,sampid, genes)
            alldats.append(res)
        newcnv = pd.concat(alldats)
        newcnv.to_csv('/tmp/pancpdo_copy_number.csv.gz',compression='gzip',index=False)
            
    if args.mutation:
        wes = sc.tableQuery("select * from syn64608378 where parentId='syn64608263'").asDataFrame()
        alldats = []
        ##go through and get every mutation file
        for index,row in wes.iterrows():
            sname = row['name'].split('--')[0]
            sid = row.id
            print(sid,sname)
            if sname in set(tsamps.other_id):
                print(sname+' in transcriptomics, using that id')
                sampid = tsamps.loc[tsamps.other_id==sname]['improve_sample_id'].values[0]
                missingsamples.append('mutation,trans,'+sname)
            elif sname in set(samps.other_id):
                print(sname+' in samples but not transcriptomics, using other id')
                sampid = samps.loc[samps.other_id==sname]['improve_sample_id'].values[0]
                missingsamples.append('mutation,notrans,'+sname)
            else:
                print('Missing sample id for '+sname)
                missingsamples.append('mutation,'+sname)
                continue
            path = sc.get(sid).path
            sampid = samps.loc[samps.other_id==sname]['improve_sample_id'].values[0]
            res = parseMutFile(path,sampid, genes)
            alldats.append(res)
        newmut = pd.concat(alldats)
        newmut.to_csv("/tmp/pancpdo_mutations.csv.gz",compression='gzip',index=False)
    pd.DataFrame(missingsamples).to_csv('missing.csv',index=False,quoting=None,header=False)
if __name__=='__main__':
    main()
