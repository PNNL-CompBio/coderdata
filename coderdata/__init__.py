from .download.downloader import download
from .dataset.dataset import (
    Dataset,
    load,
    format,
    train_test_validate,
    split_train_other,
    split_train_test_validate
)

from ._version import __version__
from ._version import __version_tuple__


from .utils.utils import version
from .utils.utils import list_datasets

try:
    import matplotlib
    import seaborn as sns
except ModuleNotFoundError:
    pass
else:
    from .utils.stats import summarize_response_metric
    from .utils.stats import plot_response_metric
    from .utils.stats import plot_2d_respones_metric