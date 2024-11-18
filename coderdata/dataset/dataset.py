from pathlib import Path
import sys

import pandas as pd

class Dataset:

    def __init__(
            self,
            name: str,
            transcriptomics: pd.DataFrame=None,
            proteomics: pd.DataFrame=None,
            mutations: pd.DataFrame=None,
            copy_number: pd.DataFrame=None,
            samples: pd.DataFrame=None,
            drugs: pd.DataFrame=None,
            mirna: pd.DataFrame=None,
            experiments: pd.DataFrame=None,
            methylation: pd.DataFrame=None,
            metabolomics: pd.DataFrame=None,
            genes: pd.DataFrame=None,
            full: pd.DataFrame=None,
            ):
        
        self.name = name
        self.transcriptomics = transcriptomics
        self.proteomics = proteomics
        self.mutations = mutations
        self.copy_number = copy_number
        self.samples = samples
        self.drugs = drugs
        self.mirna = mirna
        self.experiments = experiments
        self.methylation = methylation
        self.metabolomics = metabolomics
        self.genes = genes
        self.full = full

    
    #-----------------------------
    # getters / setters & deleters
    # ----------------------------


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
    def full(self):
        return self._full
    
    @full.setter
    def full(self, value):
        self._full = value

    @full.deleter
    def full(self):
        del self._full

    
    # ----------------------------
    # Class and Instance functions
    # ----------------------------

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


# ---------------------------------------------------------------------
# Functions that are Dataset related but not Class / Instance functions
# ---------------------------------------------------------------------

def load(name: str, directory: str|Path=Path.cwd()) -> Dataset:
    print("Processing Data...", file=sys.stderr)

    if type(directory) is not Path:
        try:
            directory = Path(directory)
            if not directory.exists():
                raise OSError(
                    f"Given path / directory does not exist: '{directory}'"
                )
        except TypeError:
            raise TypeError(
                f"Invalid path / directory defined: '{directory}'"
            )

    accepted_file_endings = ('.csv', '.tsv', '.csv.gz', '.tsv.gz')
    dataset = Dataset(name)

    for child in directory.iterdir():
        if child.name in ["genes.csv", "genes.csv.gz"]:
            dataset.genes = _load_file(child)
            print("Loaded genes dataset.", file=sys.stderr)

        if (
            child.name.startswith(name)
            and child.name.endswith(accepted_file_endings)
            ):

            dataset_type = child.name[len(name)+1:].split('.')[0]
            print(dataset_type)
            if hasattr(dataset, dataset_type):
                setattr(dataset, dataset_type, _load_file(child))

    return dataset


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


def _determine_delimiter(file_path):
    print(file_path.suffixes)
    if '.tsv' in file_path.suffixes:
        return '\t'
    else:
        return ','
