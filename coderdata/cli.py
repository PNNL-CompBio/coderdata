"""
Command Line Interface to retrieve coderdata datasets.
"""

import argparse
from os import PathLike
from pathlib import Path
from typing import Union
import sys

from .download.downloader import download as download_datasets
from .utils import version
from .utils import list_datasets

def main():
    """
    Main method containing the argument parsing and execution of
    individual subroutines.
    """
    parser = argparse.ArgumentParser(prog='coderdata', add_help=True)
    parser.set_defaults(func=info)

    # subparser for the 'download' subroutine
    subparsers = parser.add_subparsers(
        title="commands",
        dest='command'
        )
    parser_download = subparsers.add_parser(
        'download',
        help='subroutine to download datasets. See "coderdata download -h" '
             'for more options.'
        )
    parser_download.add_argument(
        '-n', '--name',
        dest='DATASET_NAME',
        type=str,
        default='all',
        help='name of the dataset to download (e.g., "beataml"). '
             'Alternatively, "all" will download the full repository of '
             'coderdata datasets. See "coderdata --list" for a complete list '
             'of available datasets. Defaults to "all"'
        )
    parser_download.add_argument(
        '-p', '--local_path',
        dest="LOCAL_PATH",
        type=check_folder,
        default=Path.cwd(),
        help='defines the folder the datasets should be stored in. Defaults '
             'to the current working directory if omitted.'
    )
    parser_download.add_argument(
        '-o', '--overwrite',
        dest="OVERWRITE",
        default=False,
        action='store_true',
        help='allow dataset files to be overwritten if they already exist.'
    )
    parser_download.set_defaults(func=download)
    
    # argument group that contains flags for additional information
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument(
        '-l', '--list',
        dest="LIST",
        action='store_true',
        help="prints list of available datasets and exits program."
    )
    grp.add_argument(
        '-v', '--version',
        dest="VERSION",
        action='store_true',
        help='prints the versions of the coderdata API and dataset and exits '
             'the program'
    )

    # checks if 'coderdata' was executed without additional arguments
    # and if so prints help message and exits
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    # parse arguments and execute defined functions by `set_default()`
    # according to which subcommands / arguments where passed on the 
    # command line
    args = parser.parse_args()
    args.func(args)


def info(args):
    """
    Helper function that takes the parsed command line arguments and 
    prints either verison information or information on the available 
    datasets depending on the arguments in ``args``.

    Parameters
    ----------
    args : Namespace
        A Namespace object that contains commandline arguments parsed by
        ``ArgumentParser.parse_args()``.
    """

    # retrieve the dataset information stored in dataset.yml via
    # coderdata.utils.list_dataset() and print the information to 
    # sys.stdout
    if args.LIST:
        print(
            '\n'
            'Available datasets\n'
            '------------------\n'
        )
        list_datasets()
        print(
            '\n'
            '------------------\n\n'
            'To download individual datasets run "coderdata download --name '
            'DATASET_NAME" where "DATASET_NAME" is for example "beataml".'
        )
    
    # retrieve version number information stored in dataset.yml via
    # coderdata.utils.version() and print the information to sys.stdout
    elif args.VERSION:
        version_numbers = version()
        print(
            *(
                f"package version: {version_numbers['package']}",
                f"dataset version: {version_numbers['dataset']}"
            ),
            sep='\n',
            file=sys.stdout,
        ) 


def download(args):
    """
    Wrapper function to download datasets via ``coderdata.download()``.
    Function passes commandline arguments to the internal download
    function.

    Parameters
    ----------
    args : Namespace
        A Namespace object that contains commandline arguments parsed by
        ``ArgumentParser.parse_args()``.
    """
    download_datasets(
        name=args.DATASET_NAME,
        local_path=args.LOCAL_PATH,
        exist_ok=args.OVERWRITE,
        )


def check_folder(path: Union[str, PathLike, Path]) -> Path:
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


# Routine to execute the main function.
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
    
    
