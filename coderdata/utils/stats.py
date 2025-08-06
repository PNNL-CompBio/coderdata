"""
Collection of helper scripts to generate general statistics on the data
contained in a CoderData Object.
"""


from copy import deepcopy

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import seaborn as sns

import coderdata as cd

def plot_2d_respones_metric(
        data: cd.Dataset,
        metric1: str,
        metric2: str,
        **kwargs: dict
    ) -> None:

    data_plot = _prepare_2d_hist_data(
        data=data.experiments,
        metrics = [metric1, metric2],
        )
    
    joint_bins = kwargs.get('joint_bins', 50)
    marginal_bins = kwargs.get('marginal_bins', 50)

    sns.jointplot(
        data=data_plot,
        x=metric2, 
        y=metric1,
        kind="hist",
        joint_kws=dict(bins=joint_bins),
        marginal_kws=dict(bins=marginal_bins)
        )

def plot_response_metric(
        data: cd.Dataset,
        metric: str='auc',
        ax: Axes=None,
        **kwargs: dict
    ) -> None:
    """
    Creates a histogram detailing the distribution of dose response 
    values for a given dose respones metric.

    If used in conjunction with `matplotlib.pyplot.subplot` or 
    `matplotlib.pyplot.subplots` and the axes object is passed to the
    function, the function populates the axes object with the generated
    plot.
    
    Parameters
    ----------
    data : coderdata.DataLoader
        A full CoderData object of a dataset
    metric : str, default='auc'
        A string that defines the response metric that should be plotted
    ax : matplotlib.axes.Axes, default=None
        An `Axes` object can be defined. This is uesful if a multipannel 
        subplot has been defined prior via `matplotlib.pyplot.subplots`.
        Passing the location of the axes to the function will then 
        populate the subplot at the given location with the generated 
        plot.
    **kwargs : dict, optional
        Additional keyword arguments that can be passed to the function
        - bins : int - sets the number of bins; passed to 
        `seaborn.histplot`
        - title : str - sets the title of the axes
        - kde : bool - adds a kernel density estimate plot into the
        histogram

    Returns
    -------
    None

    Example
    -------
    In a Jupyter Notebook environment the following snippet can be used 
    to display a histgram detailing the distribution of drug response
    AUC measures in the beataml dataset. 

    >>> import coderdata as cd
    >>> beataml = cd.DataLoader('beataml')
    >>> cd.plot_response_metric(data=beataml, metric='auc', bin=10)

    For generating multipanel plots we can make use of matplotlib and 
    the `ax` parameter of this function. Furthermore, other features / 
    parameters of the cerated figure can be changed (e.g. the title of 
    the figure via `suptitle()`). Finally it can be saved.

    >>> import coderdata as cd
    >>> import matplotlib.pyplot as plt
    >>> beataml = cd.DataLoader('beataml')
    >>> fig, axs = plt.subplots(ncols=2, figsize=(10, 5))
    >>> plot_response_metric(
    ...     data=beataml,
    ...     metric='auc', 
    ...     bins=10,
    ...     ax=axs[0]
    ...     )
    >>> plot_response_metric(
    ...     data=beataml,
    ...     metric='aac', 
    ...     bins=10,
    ...     ax=axs[0]
    ...     )
    >>> fig.set_layout_engine('tight')
    >>> fig.suptitle('Distribution of drug response values')
    >>> fig.savefig('figure.png')
    """

    # assinging values to variables based on **kwargs and defining
    # default values if not present in **kwargs
    bins_ = kwargs.get('bins', 10)
    title_ = kwargs.get('title', None)
    kde_ = kwargs.get('kde', False)

    # retrieving the data/values necessary to generate the figure
    metrics = (
        data.experiments # getting the experiments DF from the dataset
            .groupby('dose_response_metric') # grouping for later
    )
    metric_ = metrics.get_group(metric) # retrieving the desired group
    x = metric_['dose_response_value'] # getting the values

    sns.set_theme(palette='colorblind')
    p = sns.histplot(data=x, kde=kde_, bins=bins_, ax=ax)
    p.set_xlabel(metric)
    p.set_title(title_)


def summarize_response_metric(data: cd.Dataset) -> pd.DataFrame:
    """
    Helper function to extract basic statistics for the `experiments`
    object in a CoderData object. Uses `pandas.DataFrame.describe()` 
    internally to generate count, mean, standard deviation, minimum, 
    25-, 50- and 75-percentile as well as maximum for 
    `dose_response_value` for each `dose_response_metric` present in 
    `experiments`. 

    Parameters
    ----------
    data : coderdata.cd.Dataset
        A full CoderData object of a dataset

    Returns
    -------
    pandas.DataFrame
        A `pandas.DataFrame` containing basic statistics for each
        dose response metric.

    Example
    -------

    The Example assumes that a dataset with the prefix 'beataml' has 
    been downloaded previously. See also ``coderdata.download()``

    >>> import coderdata as cd
    >>> beataml = cd.DataLoader('beataml')
    >>> summary_stats = summarize_response_metric(data=beataml)
    >>> summary_stats
                            count          mean           std
    dose_response_metric                                                          
    aac                   23378.0  3.028061e-01  1.821265e-01  ...  
    auc                   23378.0  6.971939e-01  1.821265e-01  ...   
    dss                   23378.0  3.218484e-01  5.733492e-01  ...   
    ...                   ...      ...           ...           ...     
    """
    df_ret = (
        data.experiments # get experiments DF
            .groupby('dose_response_metric') # grouping by metric
            ['dose_response_value'] # value to summarize
            .describe() # get count, mean, std, etc.
            ) 

    return df_ret


def _prepare_2d_hist_data(
        data: pd.DataFrame,
        metrics: list[str]=[
            "aac", "auc", "dss",
            "fit_auc", "fit_ec50", "fit_ec50se",
            "fit_einf", "fit_hs", "fit_ic50",
            "fit_r2",
            ],
        r2: float=None,
    ) -> pd.DataFrame:
    

    metric_groups = data.groupby('dose_response_metric')
    
    if r2 is not None:
        r2_ = deepcopy(metric_groups.get_group("fit_r2"))
        r2_.rename(columns={"dose_response_value": "r2_thresh"}, inplace=True)
        r2_.drop(
            columns=[
                'source', 'time_unit', 'dose_response_metric'
                ],
            inplace=True
        )
    # print(metric_groups)
    d_ret = deepcopy(metric_groups.get_group(metrics[0]))
    d_ret.rename(columns={"dose_response_value": metrics[0]}, inplace=True)
    d_ret.drop(columns=["dose_response_metric"], inplace=True)
    
    
    for metric in metrics[1:]:
        m = deepcopy(metric_groups.get_group(metric))
        m.rename(columns={"dose_response_value": metric}, inplace=True)
        m.drop(
            columns=[
                'source', 'time_unit', 'dose_response_metric'
                ],
            inplace=True
            )

        d_ret = d_ret.merge(m, on=["improve_drug_id", "improve_sample_id", "time", "study"])

    if r2 is not None:
        d_ret = d_ret.merge(r2_, on=["improve_drug_id", "improve_sample_id", "time", "study"])
        d_ret = d_ret[d_ret["r2_thresh"] > float(r2)]
        d_ret.drop(columns=["r2_thresh"], inplace=True)


    return d_ret