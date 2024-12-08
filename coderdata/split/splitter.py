"""
Contains functions that help in generating training, testing and
validation sets for CoderData Objects.
"""

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
    training, testing and validating machine learning algorithms. 
    
    The size of the splits can be adjusted to be different from 80:10:10
    (the default)for train:test:validate. The function also allows for 
    additional optional arguments, that define the type of split that is
    performed ('mixed-set', 'drug-blind', 'cancer-blind'), if the splits 
    should be stratified (and which drug response metric to use), as
    well as a random seed to enable the creation of reproducable splits.
    Furhermore, a list of keyword arguments can be defined that will be 
    passed to the stratification function if so desired.

    Parameters
    ----------
    data : DatasetLoader
        CoderData object containing a full dataset either downloaded
        from the CoderData repository (see also 
        `coderdata.download.downloader.download_data_by_prefix`) or
        built locally via the `build_all` process. The object must first
        be loaded via `coderdata.load.loader.DatasetLoader`.
    split_type : {'mixed-set', 'drug-blind', 'cancer-blind'}, \
        default='mixed-set'

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
    ratio : tuple[int, int, int], default=(8,1,1)
        Defines the size ratio of the resulting test, train and
        validation sets. 
    stratify_by : str | None, default=None
        Defines if the training, testing and validation sets should be 
        stratified. Any value other than None indicates stratification 
        and defines which drug response value should be used as basis 
        for the stratification. _None_ indicates that no stratfication
        should be performed.
    random_state : int | RandomState | None, defaul=None
        Defines a seed value for the randomization of the splits. Will
        get passed to internal functions. Providing the seed will enable
        reproducability of the generated splits.
    **kwargs
        Additional keyword arguments that will be passed to the function
        that generates classes for the stratification 
        (see also ``_create_classes``).
    
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

    # StratifiedShuffleSplit is similar to ShuffleSplit with the added 
    # functionality to also stratify the splits according to defined 
    # class labels.
    # 
    # StratifiedShuffleSplit will be used for stratified mixed-set 
    # train/test/validate sets.

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

    # StratifiedGroupKFold generates K folds that take the group into 
    # account when generating folds, i.e. a group will only be present
    # in one fold. It further tries to stratify the folds based on the 
    # defined classes. 
    # 
    # StratifiedGroupKFold will be used for stratified drug-/sample-
    # blind splitting. 
    # 
    # The way the K folds are utilized is to combine i, j, & k folds 
    # (according to the defined ratio) into training, testing and 
    # validation sets.
    sgk = StratifiedGroupKFold(
        n_splits=sum(ratio),
        shuffle=True,
        random_state=random_state
    )

    # The "actual" splitting logic using the defined Splitters as above
    # follows here starting with the non-stratified splitting:
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
    
    # The following block contains the stratified splitting logic
    else:
        # First the classes that are needed for the stratification are
        # generated. `num_classes`, `thresh` and `quantiles` were 
        # previously defined as possible keyword arguments.
        df_full = _create_classes(
            data=df_full,
            metric=stratify_by,
            num_classes=num_classes,
            thresh=thresh,
            quantiles=quantiles,
            )
        if split_type == 'mixed-set':
            # Using StratifiedShuffleSplit to generate randomized train 
            # and 'other' set, since there is no need for grouping.
            idx_train, idx_other = next(
                sss_1.split(X=df_full, y=df_full['split_class'])
            )
            df_train = df_full.iloc[idx_train]
            df_other = df_full.iloc[idx_other]
            # Splitting 'other' further into test and validate
            idx_test, idx_val = next(
                sss_2.split(X=df_other, y=df_other['split_class'])
                )
            df_test = df_other.iloc[idx_test]
            df_val = df_other.iloc[idx_val]
        
        # using StratifiedGroupKSplit for the stratified drug-/sample-
        # blind splits.
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

            # StratifiedGroupKSplit is setup to generate K splits where
            # K=sum(ratios) (e.g. 10 if ratio=8:1:1). To obtain three 
            # sets (train/test/validate) the individual splits need to 
            # be combined (e.g. k=[1:8] -> train, k=9 -> test, k=10 -> 
            # validate). The code block below does that by combining
            # all indices (row numbers) that go into individual sets and
            # then extracting and adding those rows into the individual
            # sets.
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

    Parameters
    ----------
    data : DatasetLoader
        CoderData object containing the "full" dataset, e.g. the dataset
        that splits are based on.
    split : pandas.DataFrame
        Contains a subset of rows from the data.experiments DataFrame 
        that correspond to the generated split.

    Returns
    -------
    DatasetLoader
        A CoderData object that is a subset of ``data`` containing only
        the data points that pertain to the information in ``split``. 
    """

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
    #
    # Datapoints in the individual datasets in the CoderData object are
    # filtered by the following "rules":
    #
    # cd.drugs -> reduce based on improve_drug_id
    # cd.mutations -> reduce based on improve_sample_id
    # cd.proteomics -> reduce based on improve_sample_id
    # cd.samples -> reduce based on improve_sample_id
    # cd.transcriptomics -> reduce based on improve_sample_id
    
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

    Parameters
    ----------
    data : pandas.DataFrame
        The DataFrame containing drug response data (experiment) which 
        is subject to binning / creating classes based on the
        drug response metric
    metric : str
        The drug response metric upon which the class generation should
        be based on. Needs to be in data['drug_response_metric'].
    num_classes : int, default=2
        Number of classes that should be generated. Defaults to 2.
    quantiles : bool, default=True
        Defines whether the individual bins should be based on quantiles
        (or percentiles) instead of "evenly" spaced. If true then the
        "bin" size will be chosen such that roughly the same number of 
        data points fall into each class
    thresh : float, default=None
        Optional argument that defines a threshold other than the mean 
        of the drug response metric if ``num_classes=2``. Can be used to 
        generate "uneven" classes.

    Returns
    -------
    pandas.DataFrame
        DataFrame that is the same as the input with additional column
        that defines the established class association of each data
        point.

    Raises
    ------
    ValueError
        If the chosen ``metric`` is not present in the
        `drug_response_metric` column of ``data``.
    ValueError
        If ``num_classes`` < 2.
    ValueError
        If ``thresh`` is defined but ``num_classes`` > 2.
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
