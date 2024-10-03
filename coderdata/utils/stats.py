"""
Collection of helper scripts to generate general statistics on the data
contained in a CoderData Object.
"""

from coderdata import DatasetLoader
import pandas as pd


def summarize_response_metric(data: DatasetLoader) -> pd.DataFrame:
    
    df_ret = (
        data.experiments # get experiments DF
            .groupby('dose_response_metric') # grouping by metric
            ['dose_response_value'] # value to summarize
            .describe() # get count, mean, std, etc.
            ) 

    return df_ret

