import wget
import os

print(dir(os))

# need to downloaded Manifest file before using this script

#change directory to working directory
# os.chdir('/Users/schw939/candle_stuff/GDC')

# download the gdc client tool
def download_gdc_tool(url):
    filename = wget.download(url)
    return filename

#download_gdc_tool('https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_OSX_x64.zip')


# use the gdc client tool to download the files in your manifest
def use_gdc_tool():
    # create a new folder in the directory to store the files
        #need to make this a loop in case manifest_files is already a folder
    os.makedirs('manifest_files')

    # changing directory path so files are downloaded into proper folder
                #need to change this path so it automatically changes paths to the folder (manifest_files)
    os.chdir('./manifest_files')

    # using the gdc tool to download files
    p1 = subprocess.run(['./gdc-client', 'download', '-m', 'path/gdc_manifest.2023-07-06.txt'])
    print(p1.stdout)

#use_gdc_tool()

#access the files and clean up the dataframes to prepare for combining and mapping
def get_clean_files():
    manifest_folders = []
    sample_filenames = []
    file_paths = []
    # access the folder with all of the downloaded file folders
    for folder in os.listdir('./manifest_files'):
        manifest_folders.append(folder)
    # if this exists in the list you created, get it out
    if '.DS_Store' in manifest_folders:
        manifest_folders.remove('.DS_Store')
    # get the names of each of the files that were downloaded
    for folder_name in manifest_folders:
        for x in os.listdir('./manifest_files/' + folder_name):
            # get rid of the 'logs' folders that are contained in each folder
            if len(x) > 4:
                sample_filenames.append(x)
    # remove checkpoint string if it exists
    for each_sample in sample_filenames:
        if '.ipynb_checkpoints' in each_sample:
            sample_filenames.remove('.ipynb_checkpoints')

    # set a list that will hold all of the cleaned up dataframes
    all_dataframes = []
    # get the final file path for each file by looping through each of the lists
    for folder_n, sample in zip(manifest_folders, sample_filenames):
        # read in all of the dataframes and add them to the df list
        dataframe = pd.read_csv(
            (os.path.join('./manifest_files/' + folder_n + '/' + sample)), delimiter='\t')
        # clean up files so they read better
        dataframe = dataframe.reset_index()
        dataframe.columns = dataframe.iloc[0]
        # first rows of files don't look like they matter (tpm data is NAN)
        dataframe = dataframe[5:]
        # Grabbing the columns that we want
        if 'tpm_unstranded' in dataframe.columns:
            dataframe = dataframe[['gene_id', 'gene_name', 'tpm_unstranded']]
            all_dataframes.append(dataframe)

    return all_dataframes

#get_clean_files()

# pull/ map the info for each of these samples and then add them to one big dataset
def map_and_combine(dataframe_list):
    df_list = []
    # read in genes.csv from cell_line data to map gene_name to entrez_id's
    genes = pd.read_csv('https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/data/genes.csv')
    # read in mutations.csv from cell_line data to map improve_id to entrez_id's
    mut = pd.read_csv('https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/data/mutations.csv')

    # for each dataframe map the entrez id using the gene_name column
    for df in dataframe_list:
        mapped_entrez = (df.merge(genes, left_on='gene_name', right_on='gene_symbol', how='left').reindex(
            columns=['gene_id', 'gene_name', 'gene_type',
                     'tpm_unstranded', 'entrez_id']))
        # for each dataframe map the improve_sample_id using the entrez_id column
        mapped_improve = (mapped_entrez.merge(mut, left_on='entrez_id', right_on='entrez_id', how='left').reindex(
            columns=['gene_id', 'gene_name', 'gene_type',
                     'tpm_unstranded', 'entrez_id', 'improve_sample_id']))
        # merge all dataframes into one
        df_list.append(mapped_improve)

    final_dataframe = pd.concat(df_list)

    return final_dataframe

#map_and_combine(get_clean_files())

# getting rid of decimal places so this column can map to other id columns in datasets if needed
def clean_gene_id(final_dataframe):
    other_id = []
    for id_name in final_dataframe.gene_id:
        if '.' not in id_name:
            other_id.append(id_name)
        else:
            head, sep, tail = id_name.partition('.')
            other_id.append(head)

    final_dataframe.gene_id = other_id
    return final_dataframe

#clean_gene_id(map_and_combine(get_clean_files()))