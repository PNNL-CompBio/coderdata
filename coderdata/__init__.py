from .download.downloader import download
from .dataset.dataset import (
    Dataset,
    load,
    format,
    train_test_validate
)

# '_version.py' will be generated by hatchling once the switch away from
# setuptools.py is finished
try:
    from ._version import __version__
except ImportError:
    __version__ = '0.1.40'
try:
    from ._version import __version_tuple__
except ImportError:
    __version_tuple__ = (0, 1, 40)

from .utils.utils import version
from .utils.utils import list_datasets