import synapseclient
import pandas as pd
import argparse


def get_bladder_pdo_experiments(synObject, samples, drugs):
    # get list of syn id info
    files = list(synObject.getChildren(parent='syn64765430', includeTypes=['file']))
    # load sample sheet and format _ to .
    # load conversion table and remove trailing _T
    conversion_table = synObject.get(files[50]['id'])
    conversion_table_df = pd.read_excel(conversion_table.path, header=None)
    conversion_table_df[2] = conversion_table_df[1].str.rsplit("_", expand=True)[0]#.replace(".", "_")
    conversion_table_df[3] = conversion_table_df[2].str.replace(".", "_")
    #print(conversion_table_df.head)

    # initiate empty pd.dat.frame
    drug_df = pd.DataFrame()
    # for each drug, 
    for i in range(len(files)-4):
        drug_table_syn =synObject.get(files[i]['id'])
        drug_table = pd.read_csv(drug_table_syn.path, sep="\t")
    # melt
    # link to conversion table
    # link to sample sheet
    # Rename, add columns
        melted_single_drug = drug_table.melt(id_vars = 'Unnamed: 0', value_vars = drug_table.columns[1:], var_name="sample")
        melted_single_drug['linkID'] = melted_single_drug['sample'].str.split(".", expand=True)[0] 
        drugdata = melted_single_drug.merge(conversion_table_df, left_on = 'linkID', right_on = 0, how='left')[['Unnamed: 0', 'value', 3]]
        #print(drugdata.head)
        drugdata_with_improvesample = drugdata.merge(samples, left_on = 3, right_on='common_name')
    
    #   print(drugdata_with_improvesample.head)
        drugdata_with_improvesample = drugdata_with_improvesample[['Unnamed: 0', 'value', 'improve_sample_id']]
        #print(drugdata_with_improvesample.columns)
        drugdata_with_improvesample = drugdata_with_improvesample.rename({"Unnamed: 0" : "DOSE", 'value' : 'GROWTH'}, axis=1)
        #print(drugdata_with_improvesample.columns)

        selected_drugdata = drugdata_with_improvesample
        selected_drugdata['chem_name'] = files[i]['name'].split(")")[1].split("(")[0].split(".")[0].strip().lower()
        #print(selected_drugdata.head)
        drugdata_with_both_improveIds = selected_drugdata.merge(drugs[['improve_drug_id', 'chem_name']], how='left')
        final_drugdata = drugdata_with_both_improveIds[['DOSE', 'GROWTH', 'improve_sample_id', 'improve_drug_id']]
        final_drugdata = final_drugdata.rename({'improve_drug_id' : "Drug"}, axis=1)
        final_drugdata['study'] = "Lee etal 2018 Bladder PDOs"
        final_drugdata['source'] = "Synapse"
        final_drugdata['time'] = 6
        final_drugdata['time_unit'] = 'days'
        #print(final_drugdata.head)
        # append to dataframe
        dose_resp_df = pd.concat([drug_df, final_drugdata])
    
    return dose_resp_df


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--token', help='Synapse authentication token')
    parser.add_argument('-s', '--curSampleFile', help='Sample mapping file for bladder pdo samples')
    parser.add_argument('-d', '--drugfile', help='Drug mapping file for bladder pdo samples')
    parser.add_argument('-o', '--output', default = '/tmp/bladderpdo_doserep.tsv',help='Output file to be read into curve fitting code')

    args = parser.parse_args()
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)
    drug_df = pd.read_csv(args.drugfile, sep='\t')
    samples_df = pd.read_csv(args.curSampleFile)

    doseresponse_data = get_bladder_pdo_experiments(synObject, samples_df, drug_df)
    doseresponse_data.to_csv(args.output, sep='\t')

