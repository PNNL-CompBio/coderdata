import argparse
from .download.downloader import download_data_by_prefix

def main():
    parser = argparse.ArgumentParser(prog='coderdata')
    subparsers = parser.add_subparsers(dest='command')

    # Subcommand 'download'
    parser_download = subparsers.add_parser('download', help='Download datasets')
    parser_download.add_argument('--prefix', type=str, default=None,
                                 help='Prefix of the dataset to download (e.g., "hcmi"), "all", or leave empty for all files.')
    parser_download.set_defaults(func=download_data_by_prefix)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        # Check if 'prefix' argument is present and pass it directly to the function
        if 'prefix' in args:
            args.func(args.prefix)
        else:
            args.func()
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
    
    
