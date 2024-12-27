"""
Collection of small utility and helper functions. 
"""

from importlib import resources
import yaml

from .. import __version__
from .. import __version_tuple__


def version() -> dict:
    """
    Helper function that returns the version strings for the package and
    the dataset build.

    Returns
    -------
    dict
        Contains package and dataset build version.
    """
    return {
        'package' : __version__,
        'dataset' : f"{__version_tuple__[0]}.{__version_tuple__[1]}"
        }


def list_datasets(raw: bool=False) -> dict | None:
    """
    Hepler function that returns a list of available datasets including 
    a short description and additional information available.

    Parameters
    ----------
    raw : bool, default=False
        If set to True returns a yaml dictionary containing all 
        available datasets including additional information. If set to 
        false prints information to stdout and returns None.

    Returns
    -------
    dict | None
        Returns a dict containing the information if ``raw==True``,
        otherwise prints information to stdout and returns `None`.
    """
    with resources.open_text('coderdata', 'dataset.yml') as f:
        data_information = yaml.load(f, Loader=yaml.FullLoader)
    if raw:
        return data_information['datasets']
    else:
        datasets = data_information['datasets']
        for dataset in data_information:
            print(f'{dataset}: "{data_information[dataset]['description']}"')
