import argparse
import sys
from .download.downloader import download as download_datasets

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
    parser_download.set_defaults(func=download)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()
    args.func(args)

def download(args):
    download_datasets(name=args.DATASET_NAME)
if __name__ == '__main__':
    main()
    
    
