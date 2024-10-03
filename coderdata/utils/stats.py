"""
Collection of helper scripts to generate general statistics on the data
contained in a CoderData Object.
"""

from coderdata import DatasetLoader
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

def summarize_response_metric(data: DatasetLoader) -> pd.DataFrame:
    
    df_ret = (
        data.experiments # get experiments DF
            .groupby('dose_response_metric') # grouping by metric
            ['dose_response_value'] # value to summarize
            .describe() # get count, mean, std, etc.
            ) 

    return df_ret


def plot_response_metric(data: DatasetLoader, metric: str='auc'):

    metrics = data.experiments.groupby('dose_response_metric')
    metric_ = metrics.get_group(metric)
    x = metric_['dose_response_value']

    fig, ax = plt.subplots()

    ax.hist(x, bins=10, linewidth=0.5, edgecolor="white")
    ax.set_title(f"value distribution for {metric}")
    plt.show()
