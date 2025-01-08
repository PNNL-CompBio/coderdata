import argparse
from os import PathLike
from pathlib import Path
from typing import Union
import sys
from .download.downloader import download as download_datasets
from .utils import version
from .utils import list_datasets

def main():
    parser = argparse.ArgumentParser(prog='coderdata')
    subparsers = parser.add_subparsers(dest='command')

    # Subcommand 'download'
    parser_download = subparsers.add_parser('download', help='Download datasets')
    parser_download.add_argument(
        '-n', '--name',
        dest='DATASET_NAME',
        type=str,
        default='all',
        help='Name of the dataset to download (e.g., "beataml"). '
             'Alternatively "all" will download the full repository of '
             'coderdata datasets. See "coderdata --list" for a complete list '
             'of available datasets. Defaults to "all"'
        )
    parser_download.add_argument(
        '-p', '--local_path',
        dest="LOCAL_PATH",
        type=check_folder,
        default=Path.cwd(),
        help='Defines the folder the datasets should be stored in. Defaults '
             'to the current working directory if omitted.'
    )
    parser_download.add_argument(
        '-o', '--overwrite',
        dest="OVERWRITE",
        default=False,
        action='store_true',
        help='Allow dataset files to be overwritten if they already exist.'
    )
    parser_download.set_defaults(func=download)
    
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument(
        '-l', '--list',
        dest="LIST",
        action='store_true',
        help="Prints list of available datasets and exits program."
    )
    grp.add_argument(
        '-v', '--version',
        dest="VERSION",
        action='store_true',
        help="Prints the versions of the coderdata API and dataset and exits the program"
    )
    parser.set_defaults(func=info)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()
    args.func(args)

def info(args):
    if args.LIST:
        print(list_datasets())
    elif args.VERSION:
        print(version())


def download(args):
    download_datasets(
        name=args.DATASET_NAME,
        local_path=args.LOCAL_PATH,
        exist_ok=args.OVERWRITE,
        )


def check_folder(path: Union[str, PathLike, Path]) -> Path:
    """
    Helper function to check if a defined folder exists
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
    try:
        main()
    except KeyboardInterrupt:
        pass
    
    
