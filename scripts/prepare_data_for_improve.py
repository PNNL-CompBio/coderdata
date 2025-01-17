
import argparse
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
                

    # look up the row ids for all data items of each data source to 
    # create "<STUDY>_all.txt in /splits"


    # join the "meta data tables" like copynumber etc.


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
