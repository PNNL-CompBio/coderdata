
from copy import deepcopy
from typing import Literal

import numpy as np
from numpy.random import RandomState

import pandas as pd

from sklearn.model_selection import GroupShuffleSplit

from coderdata.load.loader import DatasetLoader

def train_test_validate(
        data: DatasetLoader,
        split_type: Literal[
            'mixed-set', 'drug-blind', 'cancer-blind', 'disjoint'
            ]='mixed-set',
        random_state: (int | RandomState | None)=None,
        ) -> tuple[DatasetLoader, DatasetLoader, DatasetLoader]:
    """
    Splits a `CoderData` object (see also
    `coderdata.load.loader.DatasetLoader`) into three subsets for 
    training, testing and validating machine learning algorithms. The 
    Size of the splits is fixed to 80:10:10 for train:test:validate. The
    function allows for additional optional arguments, that define the 
    type of split that is performed as well as a random seed to enable
    the creation of reproducable splits.

    Parameters
    ----------
    data : DatasetLoader
        CoderData object containing a full dataset either downloaded
        from the CoderData repository (see also 
        `coderdata.download.downloader.download_data_by_prefix`) or
        built locally via the `build_all` process. The object must first
        be loaded via `coderdata.load.loader.DatasetLoader`.
    split_type : Literal['mixed-set', 'drug-blind', 'cancer-blind', \
                'disjoint'], optional
        Defines the type of split that should be generated:
            
        - *mixed-set*: Splits randomly independent of drug / cancer 
            association of the samples. Individual drugs or cancer types
            can appear in all three splits
        - *drug-blind*: Splits according to drug association. Any sample
            associated with a drug will be unique to one of the splits.
            For example samples with association to drug A will only be 
            present in the train split, but never in test or validate.
        - *cancer-blind*: Splits according to cancer association. 
            Equivalent to drug-blind, except cancer types will be unique
            to splits.
        - *disjoint-set*: Splits according to both drug and cancer
            association. Makes sure that samples in train, test and
            validate are fully disjoint according to drug and cancer 
            association, i.e. both a given drug as well as cancer will 
            be unique to one of the tree splits.

        Defaults to *mixed-set*.
    
    random_state : int | RandomState | None, optional
        Defines a seed value for the randomization of the splits. Will
        get passed to internal functions. Providing the seed will enable
        reproducability of the generated splits.

        Defaults to None
    
    Returns
    -------
    Splits : tuple[DataLoader, DataLoader, DataLoader]
        Three `CoderData` objects, one each for train, test & validate.

    Raises
    -------
    ValueError : 
    If supplied `split_type` is not in the list of accepted values.

    """

    # Type checking split_type
    if split_type not in [
        'mixed-set', 'drug-blind', 'cancer-blind', 'disjoint'
        ]:
        raise ValueError(
            f"{split_type} not an excepted input for 'split_type'"
            )

    df_full = data.experiments.copy()

    if split_type in ['mixed-set', 'drug-blind', 'cancer-blind']:
        # GroupShuffleSplit is a method/class implemented by 
        # scikit-learn that enables creating splits that are grouped 
        # over values of a defined column. Specifically this will 
        # guarantee that samples/rows with a specific column value will 
        # only appear in one of the two created splits.
        #
        # n defines how often a train/test split is generated not the 
        # number of splits generated.
        #
        # For the first split we define train/"other" (other will be
        # split into train and validate later) being split 80/20
        gss = GroupShuffleSplit(
            n_splits=1,
            train_size=0.8,
            test_size=0.2,
            random_state=random_state
        )

        if split_type == 'mixed-set':
            # For the mixed-set we first need to define a combined 
            # drug and sample id, since there are multiple rows with 
            # identical sample, drug id pairs but recording differnt 
            # drug effect strength measures (e.g. ic50 and auc) 
            df_full['sample_id_drug_id'] = (
                df_full.improve_sample_id.astype(str)
                + '-'
                + df_full.improve_drug_id
            )
            # The split is grouped by the new column and genereated
            # gss.split generates a list of length n with n different 
            # train/test splits (see above). In our case we only 
            # generate one split, but we still need the 'next' to access
            # it. Row ids for items in train (idx2) and other (idx2) are
            # extracted.
            idx1, idx2 = next(
                gss.split(df_full, groups=df_full.sample_id_drug_id)
                )
        elif split_type == 'drug-blind':
            # same as above we just group only over the drug id
            idx1, idx2 = next(
                gss.split(df_full, groups=df_full.improve_drug_id)
                )
        elif split_type == 'cancer-blind':
            # same as above we just group only over the sample id
            idx1, idx2 = next(
                gss.split(df_full, groups=df_full.improve_sample_id)
                )
        else:
            raise Exception(f"Should be unreachable")

        # generate new DFs containing the subset of items extracted for
        # train and other
        df_train = df_full.iloc[idx1]
        df_other = df_full.iloc[idx2]
        
        # generate now Splitter to split other into test and validate
        gss = GroupShuffleSplit(
            n_splits=1,
            train_size=0.5,
            test_size=0.5,
            random_state=random_state
        )

        # follows same logic as previous splitting with the difference
        # that only "other" is sampled and split
        if split_type == 'mixed-set':
            idx1, idx2 = next(
                gss.split(df_other, groups=df_other.sample_id_drug_id)
                )
        elif split_type == 'drug-blind':
            idx1, idx2 = next(
                gss.split(df_other, groups=df_other.improve_drug_id)
                )
        elif split_type == 'cancer-blind':
            idx1, idx2 = next(
                gss.split(df_other, groups=df_other.improve_sample_id)
                )
        else:
            raise Exception(f"Should be unreachable")
        
        # extract itmes for test and validate from other based on the 
        # sampled indices
        df_test = df_other.iloc[idx1]
        df_val = df_other.iloc[idx2]
    elif split_type == 'disjoint':
        raise NotImplementedError('disjoint currently not implemented')
    else:
        raise Exception(
            f"`split_type` contains unexpected value '{split_type}'!"
            )
    
    # dropping the previously generated combined index for mixed-type
    # as it is no longer necessary and would otherwise appear in the 
    # finalized CoderData objects that are to be returned
    if split_type == 'mixed-set':
        df_full.drop('sample_id_drug_id', axis='columns', inplace=True)

    # generating filtered CoderData objects that contain only the 
    # respective data for each split
    data_train = _filter(data, df_train)
    data_test = _filter(data, df_test)
    data_val = _filter(data, df_val)
    
    return (data_train, data_test, data_val)


def _filter(data: DatasetLoader, split: pd.DataFrame) -> DatasetLoader:
    """
    Helper function to filter down the CoderData object(s) to create
    indipendent more concise CoderData objects for training, testing 
    and validation splits.
    """

    # cd.drugs -> reduce based on improve_drug_id
    # cd.mutations -> reduce based on improve_sample_id
    # cd.proteomics -> reduce based on improve_sample_id
    # cd.samples -> reduce based on improve_sample_id
    # cd.transcriptomics -> reduce based on improve_sample_id

    # extracting improve sample and drug ids from the provided split
    sample_ids = np.unique(split['improve_sample_id'].values)
    drug_ids = np.unique(split['improve_drug_id'].values)
    
    # creating a deep copy of the CoderData object such that any 
    # further operations on the object are not changing the original
    # object / data
    data_ret = deepcopy(data)

    # filtering each individual data type down by only the improve 
    # sample / drug ids that are present in the split (extracted above)
    data_ret.drugs = data_ret.drugs[
        data_ret.drugs['improve_drug_id'].isin(drug_ids)
        ]
    data_ret.mutations = data_ret.mutations[
        data_ret.mutations['improve_sample_id'].isin(sample_ids)
        ]
    data_ret.proteomics = data_ret.proteomics[
        data_ret.proteomics['improve_sample_id'].isin(sample_ids)
        ]
    data_ret.samples = data_ret.samples[
        data_ret.samples['improve_sample_id'].isin(sample_ids)
        ]
    data_ret.transcriptomics = data_ret.transcriptomics[
        data_ret.transcriptomics['improve_sample_id'].isin(sample_ids)
        ]
    data_ret.experiments = split
    
    return data_ret


def get_subset(df_full: pd.DataFrame, df_subset: pd.DataFrame):
    idx = df_subset.index
    return df_full.drop(idx, axis='index')