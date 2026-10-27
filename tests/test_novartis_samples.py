"""
Novartis sample sheet: one sample per model. Omics and drug responses are
joined on the model name, so a model listed twice got two sample ids with
identical data (NIBRX-2043 in 2.4).
"""
import importlib.util
import os
import sys
import types

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.modules.setdefault('synapseclient', types.ModuleType('synapseclient'))


def _load(filename):
    spec = importlib.util.spec_from_file_location(filename[:-3].replace('-', '_'),
                                                  os.path.join(ROOT, 'coderbuild', 'novartis', filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


samples = _load('01-samples-novartis.py')


class FakeSynapse:
    """Serves sample sheets as Synapse files from a folder."""
    def __init__(self, sheets):
        self.sheets = sheets

    def getChildren(self, parent, includeTypes):
        skipped = ['syn66276102', 'syn66276098', 'syn66477971']
        return [{'id': i} for i in list(self.sheets) + skipped]

    def get(self, syn_id):
        return types.SimpleNamespace(path=self.sheets[syn_id])


def test_one_sample_per_model(tmp_path):
    sheet = pd.DataFrame({
        'Sample ID': ['CPDM_0338X_PDX_PNR_01', 'CPDM_0338X_PDX_PNR_02', 'CPDM_0001X_PDX_01'],
        'MODEL_ORIGINATOR_ID': ['NIBRX-2043', 'NIBRX-2043', 'NIBRX-1000'],
        'CANCER_HISTOLOGY': ['Pancreas ductal adenocarcinoma'] * 2 + ['Breast carcinoma'],
    })
    path = tmp_path / 'sheet.csv'
    sheet.to_csv(path, index=False)
    out = samples.get_complete_novartis_sample_sheet(FakeSynapse({'syn1': str(path)}))

    assert out['common_name'].tolist() == ['NIBRX-2043', 'NIBRX-1000']
    assert out['other_id'].tolist() == ['CPDM_0338X_PDX_PNR_01', 'CPDM_0001X_PDX_01']
