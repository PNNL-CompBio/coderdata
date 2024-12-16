"""
Collection of small utility and helper functions. 
"""
from .. import __version__
from .. import __version_tuple__

def version():
    return {
        'package' : __version__,
        'dataset' : f"{__version_tuple__[0]}.{__version_tuple__[1]}"
        }


