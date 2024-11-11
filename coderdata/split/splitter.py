
from copy import deepcopy
from typing import Literal

import numpy as np
from numpy.random import RandomState

import pandas as pd

from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.model_selection import StratifiedShuffleSplit

from coderdata.load.loader import DatasetLoader

def train_test_validate(
        data: DatasetLoader,
        split_type: Literal[
            'mixed-set', 'drug-blind', 'cancer-blind'
            ]='mixed-set',
        ratio: tuple[int, int, int]=(8,1,1),
        stratify_by: (str | None)=None,
        random_state: (int | RandomState | None)=None,
        **kwargs: dict,
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

    # reading in the potential keyword arguments that will be passed to
    # _create_classes().
    thresh = kwargs.get('thresh', None)
    num_classes = kwargs.get('num_classes', 2)
    quantiles = kwargs.get('quantiles', True)

    # Type checking split_type
    if split_type not in [
        'mixed-set', 'drug-blind', 'cancer-blind'
        ]:
        raise ValueError(
            f"{split_type} not an excepted input for 'split_type'"
            )


    # A wide (pivoted) table is more easy to work with in this instance.
    # The pivot is done using all columns but the 'dose_respones_value'
    # and 'dose_respones_metric' as index. df.pivot will generate a 
    # MultiIndex which complicates things further down the line. To that
    # end 'reset_index()' is used to remove the MultiIndex
    df_full = data.experiments.copy()
    df_full = df_full.pivot(
        index = [
            'source',
            'improve_sample_id',
            'improve_drug_id',
            'study',
            'time',
            'time_unit'
            ],
        columns = 'dose_response_metric',
        values = 'dose_response_value'
    ).reset_index()

    # Defining the split sizes. 
    train_size = float(ratio[0]) / sum(ratio)
    test_val_size = float(ratio[1] + ratio[2]) / sum(ratio)
    test_size = float(ratio[1]) / (ratio[1] + ratio[2])
    validate_size = 1 - test_size

    # ShuffleSplit is a method/class implemented by scikit-learn that
    # enables creating splits where the data is shuffled and then
    # randomly distributed into train and test sets according to the 
    # defined ratio.
    # 
    # n_splits defines how often a train/test split is generated.
    # Individual splits (if more than 1 is generated) are not guaranteed
    # to be disjoint i.e. test sets from individual splits can overlap.
    # 
    # ShuffleSplit will be used for non stratified mixed-set splitting 
    # since there is no requirement for disjoint groups (i.e. drug / 
    # sample ids). 
    shs_1 = ShuffleSplit(
        n_splits=1,
        train_size=train_size,
        test_size=test_val_size,
        random_state=random_state
    )
    shs_2 = ShuffleSplit(
        n_splits=1,
        train_size=test_size,
        test_size=validate_size,
        random_state=random_state
    )

    # GroupShuffleSplit is an extension to ShuffleSplit that also
    # factors in a group that is used to generate disjoint train and 
    # test sets, e.g. in this particular case the drug or sample id to
    # generate drug-blind or sample-blind train and test sets. 
    # 
    # GroupShuffleSplit will be used for non stratified drug-/sample-
    # blind splitting, i.e. there is a requirement that instances from 
    # one group (e.g. a specific drug) are only present in the training 
    # set but not in the test set. 
    gss_1 = GroupShuffleSplit(
        n_splits=1,
        train_size=train_size,
        test_size=test_val_size,
        random_state=random_state
    )
    gss_2 = GroupShuffleSplit(
        n_splits=1,
        train_size=test_size,
        test_size=validate_size,
        random_state=random_state
    )
    sgk = StratifiedGroupKFold(
        n_splits=sum(ratio),
        shuffle=True,
        random_state=random_state
    )

    sss_1 = StratifiedShuffleSplit(
        n_splits=1,
        train_size=train_size,
        test_size=test_val_size,
        random_state=random_state
    )
    sss_2 = StratifiedShuffleSplit(
        n_splits=1,
        train_size=test_size,
        test_size=validate_size,
        random_state=random_state
    )

    if stratify_by is None:
        if split_type == 'mixed-set':
            # Using ShuffleSplit to generate randomized train and
            # 'other' set, since there is no need for grouping.
            idx1, idx2 = next(
                shs_1.split(df_full)
                )
        elif split_type == 'drug-blind':
            # Using GroupShuffleSplit to created disjoint train and
            # 'other' sets by drug id
            idx1, idx2 = next(
                gss_1.split(df_full, groups=df_full.improve_drug_id)
                )
        elif split_type == 'cancer-blind':
            # same as above we just group over the sample id
            idx1, idx2 = next(
                gss_1.split(df_full, groups=df_full.improve_sample_id)
                )
        else:
            raise Exception(f"Should be unreachable")

        # generate new DFs containing the subset of items extracted for
        # train and other
        df_train = df_full.iloc[idx1]
        df_other = df_full.iloc[idx2]

        # follows same logic as previous splitting with the difference
        # that only "other" is sampled and split
        if split_type == 'mixed-set':
            idx1, idx2 = next(
                shs_2.split(df_other, groups=None)
                )
        elif split_type == 'drug-blind':
            idx1, idx2 = next(
                gss_2.split(df_other, groups=df_other.improve_drug_id)
                )
        elif split_type == 'cancer-blind':
            idx1, idx2 = next(
                gss_2.split(df_other, groups=df_other.improve_sample_id)
                )
        else:
            raise Exception(f"Should be unreachable")
        
        # extract itmes for test and validate from other based on the 
        # sampled indices
        df_test = df_other.iloc[idx1]
        df_val = df_other.iloc[idx2]
    else:
        df_full = _create_classes(
            data=df_full,
            metric=stratify_by,
            num_classes=num_classes,
            thresh=thresh,
            quantiles=quantiles,
            )
        if split_type == 'mixed-set':
            # Using ShuffleSplit to generate randomized train and
            # 'other' set, since there is no need for grouping.
            idx_train, idx_other = next(
                sss_1.split(X=df_full, y=df_full['split_class'])
            )
            df_train = df_full.iloc[idx_train]
            df_other = df_full.iloc[idx_other]
            idx_test, idx_val = next(
                sss_2.split(X=df_other, y=df_other['split_class'])
                )
            df_test = df_other.iloc[idx_test]
            df_val = df_other.iloc[idx_val]
        elif split_type == 'drug-blind' or split_type == 'cancer-blind':
            
            if split_type == 'drug-blind':
                splitter = enumerate(
                    sgk.split(
                        X=df_full,
                        y=df_full['split_class'],
                        groups=df_full.improve_drug_id
                    )
                )
            elif split_type == 'cancer-blind':
                splitter = enumerate(
                    sgk.split(
                        X=df_full,
                        y=df_full['split_class'],
                        groups=df_full.improve_sample_id
                    )
                )
                
            idx_train = []
            idx_test = []
            idx_val = []

            for i, (idx1, idx2) in splitter:
                if i < ratio[0]:
                    idx_train.extend(idx2)
                elif i >= ratio[0] and i < (ratio[0] + ratio[1]):
                    idx_test.extend(idx2)
                elif i >= (ratio[0] + ratio[1]) and i < (ratio[0] + ratio[1] + ratio[2]):
                    idx_val.extend(idx2)
            df_full.drop(labels=['split_class'], axis=1, inplace=True)
            df_train = df_full.iloc[idx_train]
            df_test = df_full.iloc[idx_test]
            df_val = df_full.iloc[idx_val]
        else:
            raise Exception(f"Should be unreachable")


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

    # melt the wide table to store a long table in the experiments 
    # object of the returned CoderData object (to be consistent with the
    # initial orginial experiments table)
    split_long = pd.melt(
        split,
        id_vars=[
            'source',
            'improve_sample_id',
            'improve_drug_id',
            'study',
            'time',
            'time_unit',
            ],
        var_name='dose_response_metric',
        value_name='dose_response_value'
    )

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
    data_ret.experiments = split_long
    
    return data_ret


def _create_classes(
        data: pd.DataFrame,
        metric: str,
        num_classes: int=2,
        quantiles: bool=True,
        thresh: float=None,
        ) -> pd.DataFrame:
    """
    Helper function that bins experiment data into a number of defined 
    classes for use with Stratified Splits.

    """
    
    if metric not in data.columns:
        raise ValueError(
            f"Defined metric '{metric}' not present in data."
        )

    if num_classes < 2:
        raise ValueError(
            f"'num_classes must be >=2. 'num_classes' given: '{num_classes}'"
        )

    if thresh is None:
        if quantiles:
            data['split_class'] = pd.qcut(
                data[metric],
                q=num_classes,
                labels=False,
            )
        else:
            data['split_class'] = pd.cut(
                data[metric],
                bins=num_classes,
                labels=False
                )
    elif num_classes == 2:
        data['split_class'] = pd.cut(
            data[metric],
            bins=[float(min(data[metric])), thresh, float(max(data[metric]))],
            labels=False,
            include_lowest=True,
        )
    else:
        raise ValueError(
            f"'thresh' can only be defined if num_classes == 2. num_classes "
            f"passed to function: '{num_classes}'"
        )
    
    return data


def get_subset(df_full: pd.DataFrame, df_subset: pd.DataFrame):
    idx = df_subset.index
    return df_full.drop(idx, axis='index')