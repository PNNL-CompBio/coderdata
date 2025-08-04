from .utils import version
from .utils import list_datasets

try:
    import matplotlib
    import seaborn as sns
except ModuleNotFoundError:
    import warnings
    warnings.warn(
        "package was not availble. To use coderdata.utils.stats functions "
        "please make sure 'matplotlib' & 'seaborn' are available in the "
        "environment."
        )
else:
    from .stats import summarize_response_metric
    from .stats import plot_response_metric
    from .stats import plot_2d_respones_metric
