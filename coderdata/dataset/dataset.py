"""
_summary_

"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
import pickle
import sys
from typing import Literal
from typing import Optional
from typing import Union

import numpy as np
from numpy.random import RandomState
import pandas as pd
# import polars as pl

from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.model_selection import StratifiedShuffleSplit


@dataclass
class Split:
    train: Dataset
    test: Dataset
    validate: Dataset


class Dataset:

    data_format_params = {
        "samples": (
            "improve_sample_id", "cancer_type", "model_type", "common_name",
            "other_id", "other_names", "id_source", "species"
            ),
        "transcriptomics": (
            "improve_sample_id", "entrez_id", "transcriptomics"
            ),
        "proteomics": ("improve_sample_id", "entrez_id", "proteomics"),
        "mutations": ("improve_sample_id", "entrez_id", "mutation"),
        "copy_number": ("improve_sample_id", "entrez_id", "copy_number"),
        "methylation": ("improve_sample_id", "entrez_id", "methylation"),
        "experiments": (
            "improve_sample_id", "improve_drug_id", "dose_response_value"
            ),
        "drugs": ("improve_drug_id", "chem_name", "isoSMILES"),
        "genes": ("entrez_id", "gene_symbol", "other_id")
        }

    def __init__(
            self,
            name: str=None,
            transcriptomics: pd.DataFrame=None,
            proteomics: pd.DataFrame=None,
            mutations: pd.DataFrame=None,
            copy_number: pd.DataFrame=None,
            samples: pd.DataFrame=None,
            drugs: pd.DataFrame=None,
            drug_descriptors: pd.DataFrame=None,
            mirna: pd.DataFrame=None,
            experiments: pd.DataFrame=None,
            methylation: pd.DataFrame=None,
            metabolomics: pd.DataFrame=None,
            genes: pd.DataFrame=None,
            combinations: pd.DataFrame=None,
            ):
        """
        Load datasets of a specific type into predefined attributes of this class instance.

        Each attribute will be a Pandas DataFrame corresponding to a file with the dataset prefix.
        Attributes include transcriptomics, proteomics, mutations, etc.

        Parameters
        ----------
        name : str
            The name of the dataset that is stored in the object
        transcriptomics : pd.DataFrame, optional
            _description_, by default None
        proteomics : pd.DataFrame, optional
            _description_, by default None
        mutations : pd.DataFrame, optional
            _description_, by default None
        copy_number : pd.DataFrame, optional
            _description_, by default None
        samples : pd.DataFrame, optional
            _description_, by default None
        drugs : pd.DataFrame, optional
            _description_, by default None
        mirna : pd.DataFrame, optional
            _description_, by default None
        experiments : pd.DataFrame, optional
            _description_, by default None
        methylation : pd.DataFrame, optional
            _description_, by default None
        metabolomics : pd.DataFrame, optional
            _description_, by default None
        genes : pd.DataFrame, optional
            _description_, by default None
        full : pd.DataFrame, optional
            _description_, by default None
        """
        
        self.name = name
        self.transcriptomics = transcriptomics
        self.proteomics = proteomics
        self.mutations = mutations
        self.copy_number = copy_number
        self.samples = samples
        self.drugs = drugs
        self.drug_descriptors = drug_descriptors
        self.mirna = mirna
        self.experiments = experiments
        self.methylation = methylation
        self.metabolomics = metabolomics
        self.genes = genes
        self.combinations = combinations

    
    #-----------------------------
    # getters / setters & deleters
    # ----------------------------


    @property
    def data_format_params(self):
        return self._data_format_params


    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, value):
        self._name = value

    @name.deleter
    def name(self):
        del self._name


    @property
    def transcriptomics(self):
        return self._transcriptomics
    
    @transcriptomics.setter
    def transcriptomics(self, value):
        self._transcriptomics = value

    @transcriptomics.deleter
    def transcriptomics(self):
        del self._transcriptomics


    @property
    def proteomics(self):
        return self._proteomics
    
    @proteomics.setter
    def proteomics(self, value):
        self._proteomics = value

    @proteomics.deleter
    def proteomics(self):
        del self._proteomics


    @property
    def mutations(self):
        return self._mutations
    
    @mutations.setter
    def mutations(self, value):
        self._mutations = value

    @mutations.deleter
    def mutations(self):
        del self._mutations


    @property
    def copy_number(self):
        return self._copy_number
    
    @copy_number.setter
    def copy_number(self, value):
        self._copy_number = value

    @copy_number.deleter
    def copy_number(self):
        del self._copy_number


    @property
    def samples(self):
        return self._samples
    
    @samples.setter
    def samples(self, value):
        self._samples = value

    @samples.deleter
    def samples(self):
        del self._samples


    @property
    def drugs(self):
        return self._drugs
    
    @drugs.setter
    def drugs(self, value):
        self._drugs = value

    @drugs.deleter
    def drugs(self):
        del self._drugs


    @property
    def drug_descriptors(self):
        return self._drug_descriptors
    
    @drug_descriptors.setter
    def drug_descriptors(self, value):
        self._drug_descriptors = value

    @drug_descriptors.deleter
    def drug_descriptors(self):
        del self._drug_descriptors


    @property
    def mirna(self):
        return self._mirna
    
    @mirna.setter
    def mirna(self, value):
        self._mirna = value

    @mirna.deleter
    def mirna(self):
        del self._mirna


    @property
    def experiments(self):
        return self._experiments
    
    @experiments.setter
    def experiments(self, value):
        self._experiments = value

    @experiments.deleter
    def experiments(self):
        del self._experiments


    @property
    def methylation(self):
        return self._methylation
    
    @methylation.setter
    def methylation(self, value):
        self._methylation = value

    @methylation.deleter
    def methylation(self):
        del self._methylation


    @property
    def metabolomics(self):
        return self._metabolomics
    
    @metabolomics.setter
    def metabolomics(self, value):
        self._metabolomics = value

    @metabolomics.deleter
    def metabolomics(self):
        del self._metabolomics


    @property
    def genes(self):
        return self._genes
    
    @genes.setter
    def genes(self, value):
        self._genes = value

    @genes.deleter
    def genes(self):
        del self._genes


    @property
    def combinations(self):
        return self._combinations
    
    @combinations.setter
    def combinations(self, value):
        self._combinations = value

    @combinations.deleter
    def combinations(self):
        del self._combinations

    
    # ----------------------------
    # Class and Instance functions
    # ----------------------------

    def format(
            self,
            data_type: Literal[
                'transcriptomics', 'mutations', 'copy_number', 'proteomics',
                'experiments', 'combinations', 'drug_descriptor', 'drugs',
                'genes', 'samples',
                ],
            use_polars: bool=False,
            **kwargs: dict,
            ):
        return format(self, data_type=data_type, use_polars=use_polars, **kwargs)

    
    def train_test_validate(
        self,
        split_type: Literal[
            'mixed-set', 'drug-blind', 'cancer-blind'
            ]='mixed-set',
        ratio: tuple[int, int, int]=(8,1,1),
        stratify_by: Optional[str]=None,
        random_state: Optional[Union[int,RandomState]]=None,
        **kwargs: dict,
        ) -> Split:

        split = train_test_validate(
            data=self,
            split_type=split_type,
            ratio=ratio,
            stratify_by=stratify_by,
            random_state=random_state,
            **kwargs
        )

        return split


    def types(self) -> list:
        data_types = [
            'transcriptomics',
            'proteomics',
            'mutations',
            'copy_number',
            'samples',
            'drugs',
            'mirna',
            'experiments',
            'methylation',
            'metabolomics',
            'genes',
            ]
        data_types_present = []
        for data_type in data_types:
            if getattr(self, data_type) is not None:
                data_types_present.append(data_type)
        
        return data_types_present
    
    def save(self, path: Path) -> None:

        with open(path, 'wb') as f_path:
            pickle.dump(self, file=f_path)



# ---------------------------------------------------------------------
# Functions that are Dataset related but not Class / Instance functions
# ---------------------------------------------------------------------

def load(
        name: str,
        local_path: Union[str,Path]=Path.cwd(),
        from_pickle:bool=False
        ) -> Dataset:
    """
    _summary_

    Parameters
    ----------
    name : str
        _description_
    directory : str | Path, optional
        _description_, by default Path.cwd()

    Returns
    -------
    Dataset
        _description_

    Raises
    ------
    OSError
        _description_
    TypeError
        _description_
    """

    if type(local_path) is not Path:
        try:
            local_path = Path(local_path)
            if not local_path.exists():
                raise OSError(
                    f"Given path / directory does not exist: '{local_path}'"
                )
        except TypeError:
            raise TypeError(
                f"Invalid path / directory defined: '{local_path}'"
            )


    if not from_pickle:
        dataset = Dataset(name)
        accepted_file_endings = ('.csv', '.tsv', '.csv.gz', '.tsv.gz')
        print(f"Importing raw data ...", file=sys.stderr)
        for child in local_path.iterdir():
            if child.name in ["genes.csv", "genes.csv.gz"]:
                print(
                    f"Importing 'genes' from {child} ...",
                    end=' ',
                    file=sys.stderr
                    )
                dataset.genes = _load_file(child)
                print("DONE", file=sys.stderr)

            if (
                child.name.startswith(name)
                and child.name.endswith(accepted_file_endings)
                ):

                dataset_type = child.name[len(name)+1:].split('.')[0]
                print(
                    f"Importing '{dataset_type}' from {child} ...",
                    end=' ',
                    file=sys.stderr
                    )
                if hasattr(dataset, dataset_type):
                    setattr(dataset, dataset_type, _load_file(child))
                    print("DONE", file=sys.stderr)
        print(f"Importing raw data ... DONE", file=sys.stderr)
        return dataset

    else:
        accepted_file_endings = ('.pkl', '.pickle')
        for child in local_path.iterdir():
            if (
                child.name.startswith(name)
                and child.name.endswith(accepted_file_endings)
                ):
                print(f"Importing pickled data ...", end=' ', file=sys.stderr)
                with open(child, 'rb') as file:
                    dataset = pickle.load(file=file)
                print("DONE", file=sys.stderr)
                return dataset



def format(
        data: Dataset,
        data_type: Literal[
            'transcriptomics', 'mutations', 'copy_number', 'proteomics',
            'experiments', 'combinations', 'drug_descriptor', 'drugs',
            'genes', 'samples',
            ],
        use_polars: bool=False,
        **kwargs: dict,
        ):

    if data_type == "transcriptomics":
        if data.transcriptomics is None:
            raise ValueError(
                f"'{data_type}' attribute of Dataset cannot be 'None'"
                )
        ret = pd.pivot_table(
            data=data.transcriptomics,
            values='transcriptomics',
            index='entrez_id',
            columns='improve_sample_id'
            ).transpose()

    elif data_type == "mutations":
        if data.mutations is None:
            raise ValueError(
                f"'{data_type}' attribute of Dataset cannot be 'None'"
                )
        mutation_type = kwargs.get('mutation_type', None)
        if mutation_type is None:
            raise ValueError(
                "'mutation_type' must be defined if 'data_type'=='mutations'"
            )
        tmp = deepcopy(data.mutations[
            data.mutations['variant_classification'] == mutation_type
            ])
        tmp['exists'] = 1
        ret = pd.pivot_table(
            data=tmp,
            index='entrez_id',
            columns='improve_sample_id',
            values='exists',
            fill_value=0,
            ).transpose()

    elif data_type == "copy_number":
        if data.copy_number is None:
            raise ValueError(
                f"'{data_type}' attribute of Dataset cannot be 'None'"
            )
        copy_call = kwargs.get('copy_call', False)
        
        ret = pd.pivot_table(
            data=data.copy_number,
            index='entrez_id',
            columns='improve_sample_id',
            values='copy_number',
            aggfunc='mean',
            ).transpose()
        if copy_call:
            ret = ret.apply(
                pd.cut,
                bins = [0, 0.5210507, 0.7311832, 1.214125, 1.422233, 2],
                labels = ['deep del', 'het loss', 'diploid', 'gain', 'amp'],
                include_lowest=True
            )


    elif data_type == "proteomics":
        if data.proteomics is None:
            raise ValueError(
                f"'{data_type}' attribute of Dataset cannot be 'None'"
            )
        
        ret = pd.pivot_table(
            data=data.proteomics,
            values='proteomics',
            index='entrez_id',
            columns='improve_sample_id'
            ).transpose()

    elif data_type == "experiments":
        if data.experiments is None:
            raise ValueError(
                f"'{data_type}' attribute of Dataset cannot be 'None'"
            )
        shape = kwargs.get('shape', 'long')
        legal_shapes = ['long', 'wide', 'matrix']
        metrics = kwargs.get('metrics', ['fit_auc'])
        if type(metrics) is not list:
            metrics = [metrics]
        if shape not in legal_shapes:
            raise ValueError(
                f"'shape' has to be one of '{legal_shapes}'"
            )
        tmp = data.experiments[
            data.experiments['dose_response_metric'].isin(metrics)
        ]
        if shape == 'long':
            ret = tmp
        elif shape == 'wide':
            ret = tmp.pivot(
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
            ).reset_index().rename_axis(None, axis=1)
        elif shape == 'matrix':
            if len(metrics) > 1:
                raise ValueError(
                    f"'metrics' cannot contain more than 1 value if "
                    f"shape=='{shape}'"
                )
            ret = pd.pivot_table(
                data=tmp,
                values='dose_response_value',
                index='improve_drug_id',
                columns='improve_sample_id'
            )
        return ret

    elif data_type == "combinations":
        raise NotImplementedError(
            f"'data_type' {data_type} is currently not implemented"
        )

    elif data_type == "drug_descriptor":
        if data.drug_descriptors is None:
            raise ValueError(
                f"'{data_type}' attribute of Dataset cannot be 'None'"
            )    
        shape = kwargs.get('shape', 'long')
        legal_shapes = ['long', 'wide']
        if shape not in legal_shapes:
            raise ValueError(
                f"'shape' has to be one of '{legal_shapes}'"
            )
        drug_descriptor_type = kwargs.get('drug_descriptor_type', None)

        if drug_descriptor_type is None:
            tmp = data.drug_descriptors
        else:
            # TODO: potentially allow for list of columns to retain
            tmp = data.drug_descriptors[
                data.drug_descriptors[
                    'structural_descriptor'
                    ] == drug_descriptor_type
                ]
        if shape == 'long':
            ret = tmp
        else:
            ret = tmp.pivot(
                index = 'improve_drug_id',
                columns = 'structural_descriptor',
                values = 'descriptor_value'
            )


    elif data_type == "drugs":
        if data.drugs is None:
            raise ValueError(
                f"'{data_type}' attribute of Dataset cannot be 'None'"
            )
        ret = data.drugs
    
    elif data_type == "genes":
        if data.genes is None:
            raise ValueError(
                f"'{data_type}' attribute of Dataset cannot be 'None'"
            )
        ret = data.genes
    
    elif data_type == "samples":
        if data.samples is None:
            raise ValueError(
                f"'{data_type}' attribute of Dataset cannot be 'None'"
            )
        ret = data.samples
    
    else:
        raise ValueError(
            f"'data_type' contains unsupported value: '{data_type}'"
            )

    return ret

def train_test_validate(
        data: Dataset,
        split_type: Literal[
            'mixed-set', 'drug-blind', 'cancer-blind'
            ]='mixed-set',
        ratio: tuple[int, int, int]=(8,1,1),
        stratify_by: Optional[str]=None,
        random_state: Optional[Union[int,RandomState]]=None,
        **kwargs: dict,
        ) -> Split:
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
    Splits : Split
        A ``Split`` object that contains three Dataset objects as 
        attributes (``Split.train``, ``Split.test``,
        ``Split.validate``) 

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
            df_train = df_train.drop(labels=['split_class'], axis=1)
            df_other = df_full.iloc[idx_other]
            # Splitting 'other' further into test and validate
            idx_test, idx_val = next(
                sss_2.split(X=df_other, y=df_other['split_class'])
                )
            df_test = df_other.iloc[idx_test]
            df_test = df_test.drop(labels=['split_class'], axis=1)
            df_val = df_other.iloc[idx_val]
            df_val = df_val.drop(labels=['split_class'], axis=1)
        
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
                elif (
                        i >= (ratio[0] + ratio[1])
                        and i < (ratio[0] + ratio[1] + ratio[2])
                     ):
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
    
    return Split(data_train, data_test, data_val)



def _load_file(file_path: Path) -> pd.DataFrame:
    if file_path.suffix == '.gz':
        return pd.read_csv(
            file_path,
            compression='gzip',
            delimiter=_determine_delimiter(file_path),
            )
    else:
        return pd.read_csv(
            file_path,
            delimiter=_determine_delimiter(file_path),
            )


def _determine_delimiter(file_path: Path) -> str:
    if '.tsv' in file_path.suffixes:
        return '\t'
    else:
        return ','


def _filter(data: Dataset, split: pd.DataFrame) -> Dataset:
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
