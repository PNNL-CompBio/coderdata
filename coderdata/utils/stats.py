"""
Collection of helper scripts to generate general statistics on the data
contained in a CoderData Object.
"""

from coderdata import DatasetLoader
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

def summarize_response_metric(data: DatasetLoader) -> pd.DataFrame:
    """
    Helper function to extract basic statistics for the `experiments`
    object in a CoderData object. Uses `pandas.DataFrame.describe()` 
    internally to generate count, mean, standard deviation, minimum, 
    25-, 50- and 75-percentile as well as maximum for 
    `dose_response_value` for each `dose_response_metric` present in 
    `experiments`. 

    Parameters
    ----------
    data : coderdata.DatasetLoader
        A full CoderData object of a dataset

    Returns
    -------
    pandas.DataFrame
        A `pandas.DataFrame` containing basic statistics for each
        dose response metric.
    """
    df_ret = (
        data.experiments # get experiments DF
            .groupby('dose_response_metric') # grouping by metric
            ['dose_response_value'] # value to summarize
            .describe() # get count, mean, std, etc.
            ) 

    return df_ret


def plot_response_metric(
        data: DatasetLoader,
        metric: str='auc',
        **kwargs: dict
    ) -> Figure:
    """
    Will cerate a `matplot.figure.Figure` object and return it, e.g. to 
    be saved locally. If the only purpose is to display the figure in a
    Jupyter notebook for example the return value does not need to be
    caught.

    Parameters
    ----------
    data : coderdata.DataLoader
        A full CoderData object of a dataset
    metric : str, default='auc'
        A string that defines the response metric that should be plotted
    **kwargs : dict, optional
        Additional keyword arguments that can be passed to the function

    Returns
    -------
    matplotlib.figure.Figure
        A `Figure` object that contains the generated figure. Can be 
        passed to a variable and saved locally via `Figure.savefig()`
    """

    # assinging values to variables based on **kwargs and defining
    # default values if not present in **kwargs
    bins = kwargs.get('bins', 10)
    title = kwargs.get('title', f"Value distributino for metric '{metric}'")

    # retrieving the data/values necessary to generate the figure
    metrics = (
        data.experiments # getting the experiments DF from the dataset
            .groupby('dose_response_metric') # grouping for later
    )
    metric_ = metrics.get_group(metric) # retrieving the desired group
    x = metric_['dose_response_value'] # getting the values

    fig, ax = plt.subplots() # generating the "Plot objects"

    ax.hist(x, bins=bins, linewidth=0.5, edgecolor="white")
    ax.set_title(title)
    
    return fig
