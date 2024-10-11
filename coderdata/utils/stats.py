"""
Collection of helper scripts to generate general statistics on the data
contained in a CoderData Object.
"""


from copy import deepcopy

import numpy as np

from coderdata import DatasetLoader
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import seaborn as sns


def split_experiments_by_study(data: DatasetLoader) -> dict:
    """
    Splits the CoderData object into multiple smaller CoderData objects
    according to the `study` recorded in the ``.experiments`` table in 
    the CoderData object.

    Parameters
    ----------
    data : DatasetLoader
        The CoderData object containing the data set loaded into memory
        via ``coderdata.DatasetLoader()``.

    Returns
    -------
    dict
        A dictionary dict[study, data] where keys `study` are the names 
        of the study in the ``.experiments`` part of the imported 
        CoderData object and values `data` are the filtered smaller
        CoderData objects containing only data corresponding to the 
        study. 
    """

    df_ret = {}
    experiments = data.experiments
    
    # creating the groups based on 'study' to itterate over 
    groups = experiments.groupby('study')
    for name, group in groups:

        # extracting improve sample and drug ids from the provided split
        sample_ids = list(np.unique(group['improve_sample_id'].values))
        drug_ids = list(np.unique(group['improve_drug_id'].values))
        
        # creating new CoderData objects that contain only data
        # pertaining to the study defined by the previous grouping
        df_ret[name] = _filter(
            data=data, sample_ids=sample_ids, drug_ids=drug_ids, study=name
            )
    
    return df_ret


def _filter(
        data: DatasetLoader,
        sample_ids: list,
        drug_ids: list,
        study: str=None,
        ) -> DatasetLoader:
    """
    Helper function to filter down the CoderData object(s) to create
    independent more concise CoderData objects for further processing.
    This can be either splitting a dataset according to the different 
    drug response studies (e.g. the broad_sanger dataset) or if small 
    subsets need to be extracted (e.g. training / testing splits for 
    machine learning)

    Parameters
    ----------
    data : DatasetLoader
        Contains a full CoderData object imported/loaded via 
        ``cd.DataLoader``
    sample_ids : list
        A list of improve_sample_id[s] that the CoderData object should
        be filtered to
    drug_ids : list
        A list of improve_drug_id[s] that the CoderData object should 
        be filtered to
    study : str, default = None
        The drug response study that the CoderData object should be 
        filtered to. This argument is only important for filtering the
        broad_sanger dataset if the splitting / filtering of the data 
        set is based on the drug response study

    Returns
    -------
    DatasetLoader
        The filtered CoderData object
    
    Notes
    -----

    Different data types of the CoderData object are going to be 
    filtered using either the improve_sample_id or the improve_drug_id.
    
    - cd.copynumber -> reduce based on ``improve_sample_id``
    - cd.drugs -> reduce based on ``improve_drug_id``
    - cd.experiments -> reduce based on ``study`` (only applicable if 
      the dataset is broad_sanger)
    - cd.mutations -> reduce based on ``improve_sample_id``
    - cd.proteomics -> reduce based on ``improve_sample_id``
    - cd.samples -> reduce based on ``improve_sample_id``
    - cd.transcriptomics -> reduce based on ``improve_sample_id``
    
    """

    # creating a deep copy of the CoderData object such that any 
    # further operations on the object are not changing the original
    # object / data
    data_ret = deepcopy(data)

    # filtering each individual data type down by only the improve 
    # sample / drug ids that are present in the study
    if not data_ret.copy_number.empty:
        data_ret.copy_number = data_ret.copy_number[
            data_ret.copy_number['improve_sample_id'].isin(sample_ids)
        ]
    if not data_ret.drugs.empty:
        data_ret.drugs = data_ret.drugs[
            data_ret.drugs['improve_drug_id'].isin(drug_ids)
            ]
    if not data_ret.mutations.empty:
        data_ret.mutations = data_ret.mutations[
            data_ret.mutations['improve_sample_id'].isin(sample_ids)
            ]
    if not data_ret.proteomics.empty:
        data_ret.proteomics = data_ret.proteomics[
            data_ret.proteomics['improve_sample_id'].isin(sample_ids)
            ]
    if not data_ret.samples.empty:
        data_ret.samples = data_ret.samples[
            data_ret.samples['improve_sample_id'].isin(sample_ids)
            ]
    if not data_ret.transcriptomics.empty:
        data_ret.transcriptomics = data_ret.transcriptomics[
            data_ret.transcriptomics['improve_sample_id'].isin(sample_ids)
            ]
    if not data_ret.experiments.empty:
        data_ret.experiments = data_ret.experiments[
            data_ret.experiments['study'] == study
        ]
    # TODO: do we also need to split the gene table?
    
    return data_ret

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


def plot_response_metric(
        data: DatasetLoader,
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
