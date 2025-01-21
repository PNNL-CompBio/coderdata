
import argparse
import functools as ft
from os import PathLike
from pathlib import Path
from pathlib import PurePath
from typing import Union
import sys

import coderdata as cd
import pandas as pd

def main():

    main_parser = argparse.ArgumentParser(add_help=True)

    command_parsers = main_parser.add_subparsers(
        dest="command",
        title="commands",
        required=True,
    )

    p_shared_args = argparse.ArgumentParser(add_help=False)
    p_shared_args.add_argument(
        '-w', '--work_dir',
        dest='WORKDIR',
        type=_check_folder,
        default=Path.cwd(),
        )
    p_shared_args.add_argument(
        '--overwrite',
        dest='OVERWRITE',
        action='store_true',
        )

    p_setup_workflow = command_parsers.add_parser(
        "setup",
        parents=[p_shared_args],
        add_help=True,
    )
    p_setup_workflow.set_defaults(func=setup_workflow)

    p_download_datasets = command_parsers.add_parser(
        "download",
        parents=[p_shared_args],
        add_help=True
    )
    p_download_datasets.set_defaults(func=download_datasets)
    
    p_process_datasets = command_parsers.add_parser(
        "process",
        parents=[p_shared_args],
        add_help=True
    )
    p_process_datasets.set_defaults(func=process_datasets)
    p_process_datasets.add_argument(
        '-s', '--split_type', dest="SPLIT_TYPE",
        type=str,
        choices=['mixed-set', 'drug-blind', 'cancer-blind'],
        default='mixed-set',
    )
    p_process_datasets.add_argument(
        '-n', '--num_splits', dest='NUM_SPLITS',
        type=int,
        default=10
    )
    p_process_datasets.add_argument(
        '-g', '--gene_table', dest='GENE_TABLE',
        type=str,
        required=True
    )

    p_all = command_parsers.add_parser(
        "all",
        parents=[p_shared_args],
        add_help=True
    )
    p_all.set_defaults(func=full_workflow)

    if len(sys.argv) == 1:
        main_parser.print_help(sys.stderr)
        sys.exit(0)
    try:
        args = main_parser.parse_args()
    except FileNotFoundError as e:
        sys.exit(e)
    except ValueError as e:
        sys.exit(e)
    args.func(args)


def full_workflow(args):
    setup_workflow(args)
    download_datasets(args)


def process_datasets(args):
    
    
    local_path = args.WORKDIR.joinpath('data_in_tmp')
    
    # getting the info which datasets are available
    data_sets_info = cd.list_datasets(raw=True)
    
    # loading all available datasets into a dict where the dataset name
    # is the key
    data_sets = {}
    for data_set in data_sets_info.keys():
        data_sets[data_set] = cd.load(name=data_set, local_path=local_path)


    #-------------------------------------------------------------------
    # concatting all experiments / responses to create response.tsv
    #-------------------------------------------------------------------
    experiments = []
    for data_set in data_sets_info.keys():
        # not all Datasets have experiments / drug response data
        if data_sets[data_set].experiments is not None:
            # formatting existing response data to wide
            experiment = data_sets[data_set].format(
                data_type='experiments',
                shape='wide',
                metrics=[
                    'fit_auc',
                    'fit_ic50',
                    'fit_r2',
                    'fit_ec50se',
                    'fit_einf',
                    'fit_hs',
                    'aac',
                    'auc',
                    'dss',
                ],
            )
            experiments.append(experiment)
    
    # concatenating existing response data and "clean up"
    response_data = pd.concat(experiments, axis=0, ignore_index=True)
    # TODO: potentially more columns must be renamed
    # (e.g. fit_auc to auc). If so this would happen here
    response_data.rename(
        columns={'improve_drug_id': 'improve_chem_id'},
        inplace=True,
        )
    # exporting the drug response data to 'y_data/response.tsv'
    outfile_path = args.WORKDIR.joinpath("data_out", "y_data", "response.tsv")
    response_data.to_csv(
        path_or_buf=outfile_path,
        index=False,
        sep='\t',
        )
    # temporary addition of "index column" to serve as a reference for
    # the extraction of split files
    response_data['index'] = response_data.index

    #-------------------------------------------------------------------
    # creation of splits
    #-------------------------------------------------------------------

    splits_folder = args.WORKDIR.joinpath('data_out', 'splits')
    split_type = args.SPLIT_TYPE
    # TODO: potentially change vars to be read from `args`
    ratio = (8,1,1)
    stratify_by = None
    random_state = None

    for data_set in data_sets_info.keys():
        if data_sets[data_set].experiments is not None:

            # getting "<DATASET>_all.txt"
            drug_response_rows = (
                data_sets['mpnst']
                .experiments[
                    ['improve_sample_id', 'improve_drug_id', "time", "study"]
                    ]
                .drop_duplicates()
                )
            drug_response_rows.rename(
                    columns={'improve_drug_id': 'improve_chem_id'},
                    inplace=True,
                    )
            row_nums = pd.merge(
                response_data,
                drug_response_rows,
                how='inner',
                on=['improve_sample_id', 'improve_chem_id', "time", "study"]
                )
            outfile_path = splits_folder.joinpath(f"{data_set}_all.txt")
            row_nums.to_csv(
                path_or_buf=outfile_path,
                columns=['index'],
                index=False,
                header=False
                )

            # generation of the actual splits

            # TODO: potentially clean this up with a function that is

            splits = {}
            for i in range(0, args.NUM_SPLITS):
                splits[i] = data_sets[data_set].train_test_validate(
                    split_type=split_type,
                    ratio=ratio,
                    stratify_by=stratify_by,
                    random_state=random_state
                    )
                train_keys = (
                    splits[i]
                    .train
                    .experiments[[
                        'improve_sample_id',
                        'improve_drug_id',
                        "time",
                        "study"
                        ]]
                    .drop_duplicates()
                )
                train_keys.rename(
                    columns={'improve_drug_id': 'improve_chem_id'},
                    inplace=True,
                )
                row_nums = pd.merge(
                    response_data,
                    train_keys,
                    how='inner',
                    on=['improve_sample_id', 'improve_chem_id', "time", "study"],
                    )
                outfile_path = splits_folder.joinpath(f"{data_set}_split_{i}_train.txt")
                row_nums.to_csv(
                    path_or_buf=outfile_path,
                    columns=['index'],
                    index=False,
                    header=False
                    )
                
                test_keys = (
                    splits[i]
                    .test
                    .experiments[[
                        'improve_sample_id',
                        'improve_drug_id',
                        "time",
                        "study"
                        ]]
                    .drop_duplicates()
                )
                test_keys.rename(
                    columns={'improve_drug_id': 'improve_chem_id'},
                    inplace=True,
                )
                row_nums = pd.merge(
                    response_data,
                    test_keys,
                    how='inner',
                    on=['improve_sample_id', 'improve_chem_id', "time", "study"],
                    )
                outfile_path = splits_folder.joinpath(f"{data_set}_split_{i}_test.txt")
                row_nums.to_csv(
                    path_or_buf=outfile_path,
                    columns=['index'],
                    index=False,
                    header=False
                    )
                
                val_keys = (
                    splits[i]
                    .validate
                    .experiments[[
                        'improve_sample_id',
                        'improve_drug_id',
                        "time",
                        "study"
                        ]]
                    .drop_duplicates()
                )
                val_keys.rename(
                    columns={'improve_drug_id': 'improve_chem_id'},
                    inplace=True,
                )
                row_nums = pd.merge(
                    response_data,
                    val_keys,
                    how='inner',
                    on=['improve_sample_id', 'improve_chem_id', "time", "study"],
                    )
                outfile_path = splits_folder.joinpath(f"{data_set}_split_{i}_val.txt")
                row_nums.to_csv(
                    path_or_buf=outfile_path,
                    columns=['index'],
                    index=False,
                    header=False
                    )

    #-------------------------------------------------------------------
    # getting common / reference gene symbols
    #-------------------------------------------------------------------

    # TODO: potentially add mapping to the genes table in coderdata
    # currently we do not make use of the 'genes' DataFrame in a Dataset
    # object. The gene symbol information comes directly from HGNC.
    # There are instances where the entrez_id that is recoreded in the
    # expression / transcriptome is not in HGNC. Those currenly result
    # in NaNs for the gene symbol

    data_gene_names = pd.read_table(
        filepath_or_buffer=args.GENE_TABLE,
        )
    data_gene_names.rename(
            columns={
                'NCBI Gene ID': 'entrez_id',
                'Ensembl gene ID': 'ensemble_gene_id',
                'Approved symbol': 'gene_symbol'
                },
            inplace=True,
            )
    data_gene_names.dropna(axis=0, subset='entrez_id', inplace=True)
    data_gene_names['entrez_id'] = data_gene_names['entrez_id'].astype(int)      

    #-------------------------------------------------------------------
    # create gene expression master table
    #-------------------------------------------------------------------

    merged_transcriptomics = merge_master_tables(
        args=args,
        data_sets=data_sets,
        data_type='transcriptomics'
        )
    
    # TODO: Potentially cast 'NaN's to 0

    # merging ensemble gene id & gene symbol into the transcriptomics 
    # data
    merged_transcriptomics = pd.merge(
        merged_transcriptomics,
        data_gene_names[[
            'entrez_id',
            'ensemble_gene_id',
            'gene_symbol'
        ]],
        how='left',
        on='entrez_id',
    )

    # moving ensemble_id & gene_symbol columns to the front of the table
    # such that when transposing the DataFrame they are row 3 and 2
    # respectively
    merged_transcriptomics.insert(
        1,
        'ensemble_gene_id',
        merged_transcriptomics.pop('ensemble_gene_id')
    )
    merged_transcriptomics.insert(
        1,
        'gene_symbol',
        merged_transcriptomics.pop('gene_symbol')
    )

    # writing the expression datatable to '/x_data/*_expression.tsv'
    outfile_path = args.WORKDIR.joinpath(
        "data_out",
        "x_data",
        "cancer_gene_expression.tsv"
    )
    merged_transcriptomics.transpose().to_csv(
        path_or_buf=outfile_path,
        sep='\t',
        header=False
    )


    #-------------------------------------------------------------------
    # create copynumber master table
    #-------------------------------------------------------------------

    merged_copy_number = merge_master_tables(args, data_sets=data_sets, data_type='copy_number')

    merged_copy_number = pd.merge(
        merged_copy_number,
        data_gene_names[[
            'entrez_id',
            'ensemble_gene_id',
            'gene_symbol'
        ]],
        how='left',
        on='entrez_id',
    )

    merged_copy_number.insert(
        1,
        'ensemble_gene_id',
        merged_copy_number.pop('ensemble_gene_id')
    )
    merged_copy_number.insert(
        1,
        'gene_symbol',
        merged_copy_number.pop('gene_symbol')
    )

    # writing the expression datatable to '/x_data/*_copy_number.tsv'
    outfile_path = args.WORKDIR.joinpath(
        "data_out",
        "x_data",
        "cancer_copy_number.tsv"
    )
    merged_copy_number.transpose().to_csv(
        path_or_buf=outfile_path,
        sep='\t',
        header=False
    )
    # join the "meta data tables" like copynumber etc.


def merge_master_tables(args, data_sets, data_type: str='transcriptomics'):
    """
    Helper function to merge several DataTables into one master table

    Parameters
    ----------
    args : _type_
        _description_
    data_type : str, optional
        _description_, by default 'transcriptomics'

    Returns
    -------
    _type_
        _description_
    """

    # creating a list that contains all DataFrames to be merged
    dfs_to_merge = []
    for data_set in data_sets:
        if data_sets[data_set].experiments is not None:
            if (
                data_type in ['transcriptomics', 'copy_number'] and 
                getattr(data_sets[data_set], data_type, None) is not None
            ):
                dfs_to_merge.append(
                    data_sets[data_set].format(data_type=data_type)
                    )

    merged_data = ft.reduce(
        lambda left_df, right_df: pd.merge(
            left_df,
            right_df,
            on='entrez_id',
            how='outer',
        ),
        dfs_to_merge,
    )

    # temporary fix to values that should be int but currently aren't 
    # in the coderdata dataset storage
    if not merged_data.index.dtype == int:
        merged_data.index = merged_data.index.astype(int)

    return merged_data


def download_datasets(args):
    local_path = args.WORKDIR.joinpath('data_in_tmp')
    exist_ok = args.OVERWRITE
    try:
        cd.download(name='all', local_path=local_path, exist_ok=exist_ok)
    except FileExistsError:
        sys.exit("data files already exist")


def setup_workflow(args):
    
    # Create the folder structure according to the IMPROVE pipeline
    # including some temporary working directories. The structure will
    # look like this:
    # 
    # .
    # ├── data_in_tmp     <- will contain the downloaded datasets etc.
    # └── data_out        <- prepared data for IMPROVE pipeline
    #     ├── splits      <- contains n split files per dataset
    #     ├── x_data      <- contains combined "master tables" of data
    #     └── y_data      <- contains drug responses

    parent = args.WORKDIR
    exist_ok = args.OVERWRITE

    data_in = parent.joinpath('data_in_tmp')
    data_out = parent.joinpath('data_out')
    splits = data_out.joinpath('splits')
    x_data = data_out.joinpath('x_data')
    y_data = data_out.joinpath('y_data')

    try:
        data_in.mkdir(exist_ok=exist_ok)
        data_out.mkdir(exist_ok=exist_ok)
        splits.mkdir(exist_ok=exist_ok)
        x_data.mkdir(exist_ok=exist_ok)
        y_data.mkdir(exist_ok=exist_ok)
    except FileExistsError:
        sys.exit(
            "Some folders already exist. To ovewrite contents use "
            "commandline argument '--overwrite'"
            )


def _check_folder(path: Union[str, PathLike, Path]) -> Path:
    """
    Helper function to check if a defined folder exists.

    Returns
    -------
    Path
        Cleaned path object with the absolute path to the folder passed
        to the function.  

    Raises
    ------
    TypeError
        If passed path argument is not of the requested type.
    OSError
        If the passed path argument does not link to a valid existing
        folder.
    """

    if not isinstance(path, (str, PathLike, Path)): 
        raise TypeError(
            f"'path' must be of type str, PathLike or Path. Supplied argument "
            f"is of type {type(path)}."
        )
    if not isinstance(path, Path):
        abs_path = Path(path).absolute()
    else:
        abs_path = path.absolute()
    
    if not abs_path.is_dir():
        raise OSError(
            f"The defined folder path '{path}' does not exist or is not a "
            f"folder."
            )
    
    return abs_path

if __name__ == '__main__':
    try: main()
    except KeyboardInterrupt: pass
