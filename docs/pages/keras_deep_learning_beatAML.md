---
layout: default
title: CoderData
---

<link rel="stylesheet" href="assets/css/style.css">




```python
import pandas as pd
import numpy as np
import coderdata as cd
```


```python
import keras
from numpy import loadtxt
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from rdkit import Chem
from rdkit.Chem import AllChem
from keras.models import Model
from keras.layers import Input, Dense, concatenate
from keras.optimizers import Adam
from keras.losses import MeanSquaredError
from keras.metrics import MeanAbsoluteError
from keras import layers
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
```


```python
beataml = cd.DatasetLoader("beataml")

```

    Processing Data...



```python
beataml.info()
```

    Dataset Type: beataml
    Beat acute myeloid leukemia (BeatAML) data was collected though GitHub and Synapse.
    
    Available Datatypes and Their Formats:
    - drugs: Format not specified
    - experiments: Format not specified
    - mutations: long format
    - proteomics: long format
    - samples: long format
    - transcriptomics: long format



```python
beataml.drugs
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>improve_drug_id</th>
      <th>chem_name</th>
      <th>formula</th>
      <th>weight</th>
      <th>inCHIKey</th>
      <th>canSMILES</th>
      <th>isoSMILES</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>SMI_43</td>
      <td>17-AAG</td>
      <td>C31H43N3O8</td>
      <td>585.70</td>
      <td>AYUNIORJHRXIBJ-TXHRRWQRSA-N</td>
      <td>CC1CC(C(C(C=C(C(C(C=CC=C(C(=O)NC2=CC(=O)C(=C(C...</td>
      <td>C[C@H]1C[C@@H]([C@@H]([C@H](/C=C(/[C@@H]([C@H]...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>SMI_3222</td>
      <td>A-674563</td>
      <td>C22H22N4O</td>
      <td>358.40</td>
      <td>BPNUQXPIQBZCMR-IBGZPJMESA-N</td>
      <td>CC1=C2C=C(C=CC2=NN1)C3=CC(=CN=C3)OCC(CC4=CC=CC...</td>
      <td>CC1=C2C=C(C=CC2=NN1)C3=CC(=CN=C3)OC[C@H](CC4=C...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>SMI_146</td>
      <td>ABT-737</td>
      <td>C42H45ClN6O5S2</td>
      <td>813.40</td>
      <td>HPLNQCPCUACXLM-PGUFJCEWSA-N</td>
      <td>CN(C)CCC(CSC1=CC=CC=C1)NC2=C(C=C(C=C2)S(=O)(=O...</td>
      <td>CN(C)CC[C@H](CSC1=CC=CC=C1)NC2=C(C=C(C=C2)S(=O...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>SMI_57721</td>
      <td>AKT Inhibitor IV</td>
      <td>C31H27IN4S</td>
      <td>614.50</td>
      <td>NAYRELMNTQSBIN-UHFFFAOYSA-M</td>
      <td>CC[N+]1=C(N(C2=C1C=C(C=C2)C3=NC4=CC=CC=C4S3)C5...</td>
      <td>CC[N+]1=C(N(C2=C1C=C(C=C2)C3=NC4=CC=CC=C4S3)C5...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>SMI_57735</td>
      <td>AMPK Inhibitor</td>
      <td>C24H25N5O</td>
      <td>399.50</td>
      <td>XHBVYDAKJHETMP-UHFFFAOYSA-N</td>
      <td>C1CCN(CC1)CCOC2=CC=C(C=C2)C3=CN4C(=C(C=N4)C5=C...</td>
      <td>C1CCN(CC1)CCOC2=CC=C(C=C2)C3=CN4C(=C(C=N4)C5=C...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>35420</th>
      <td>SMI_414</td>
      <td>Venetoclax</td>
      <td>C45H50ClN7O7S</td>
      <td>868.40</td>
      <td>LQBVNQSMGBZMKD-UHFFFAOYSA-N</td>
      <td>CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC...</td>
      <td>CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC...</td>
    </tr>
    <tr>
      <th>35421</th>
      <td>SMI_1062</td>
      <td>Vismodegib</td>
      <td>C19H14Cl2N2O3S</td>
      <td>421.30</td>
      <td>BPQMGSKTAYIVFO-UHFFFAOYSA-N</td>
      <td>CS(=O)(=O)C1=CC(=C(C=C1)C(=O)NC2=CC(=C(C=C2)Cl...</td>
      <td>CS(=O)(=O)C1=CC(=C(C=C1)C(=O)NC2=CC(=C(C=C2)Cl...</td>
    </tr>
    <tr>
      <th>35422</th>
      <td>SMI_3146</td>
      <td>Volasertib</td>
      <td>C34H50N8O3</td>
      <td>618.80</td>
      <td>SXNJFOWDRLKDSF-XKHVUIRMSA-N</td>
      <td>CCC1C(=O)N(C2=CN=C(N=C2N1C(C)C)NC3=C(C=C(C=C3)...</td>
      <td>CC[C@@H]1C(=O)N(C2=CN=C(N=C2N1C(C)C)NC3=C(C=C(...</td>
    </tr>
    <tr>
      <th>35423</th>
      <td>SMI_1137</td>
      <td>XAV-939</td>
      <td>C14H11F3N2OS</td>
      <td>312.31</td>
      <td>KLGQSVMIPOVQAX-UHFFFAOYSA-N</td>
      <td>C1CSCC2=C1N=C(NC2=O)C3=CC=C(C=C3)C(F)(F)F</td>
      <td>C1CSCC2=C1N=C(NC2=O)C3=CC=C(C=C3)C(F)(F)F</td>
    </tr>
    <tr>
      <th>35424</th>
      <td>SMI_181</td>
      <td>YM-155</td>
      <td>C20H19BrN4O3</td>
      <td>443.30</td>
      <td>QBIYUDDJPRGKNJ-UHFFFAOYSA-M</td>
      <td>CC1=[N+](C2=C(N1CCOC)C(=O)C3=CC=CC=C3C2=O)CC4=...</td>
      <td>CC1=[N+](C2=C(N1CCOC)C(=O)C3=CC=CC=C3C2=O)CC4=...</td>
    </tr>
  </tbody>
</table>
<p>35425 rows × 7 columns</p>
</div>




```python
beataml.experiments

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>source</th>
      <th>improve_sample_id</th>
      <th>improve_drug_id</th>
      <th>study</th>
      <th>auc</th>
      <th>ic50</th>
      <th>ec50</th>
      <th>ec50se</th>
      <th>r2fit</th>
      <th>einf</th>
      <th>hs</th>
      <th>aac1</th>
      <th>auc1</th>
      <th>dss1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>synapse</td>
      <td>4634</td>
      <td>SMI_43</td>
      <td>BeatAML</td>
      <td>0.5716</td>
      <td>NaN</td>
      <td>11.430</td>
      <td>0.000</td>
      <td>0.0000</td>
      <td>0.5716</td>
      <td>3.525</td>
      <td>0.4284</td>
      <td>0.5716</td>
      <td>0.3649</td>
    </tr>
    <tr>
      <th>1</th>
      <td>synapse</td>
      <td>4634</td>
      <td>SMI_3222</td>
      <td>BeatAML</td>
      <td>0.6310</td>
      <td>NaN</td>
      <td>10.890</td>
      <td>0.000</td>
      <td>-0.0000</td>
      <td>0.6310</td>
      <td>3.171</td>
      <td>0.3690</td>
      <td>0.6310</td>
      <td>0.2989</td>
    </tr>
    <tr>
      <th>2</th>
      <td>synapse</td>
      <td>4634</td>
      <td>SMI_146</td>
      <td>BeatAML</td>
      <td>0.7496</td>
      <td>NaN</td>
      <td>11.590</td>
      <td>0.000</td>
      <td>0.0000</td>
      <td>0.7496</td>
      <td>3.594</td>
      <td>0.2504</td>
      <td>0.7496</td>
      <td>0.1671</td>
    </tr>
    <tr>
      <th>3</th>
      <td>synapse</td>
      <td>4634</td>
      <td>SMI_57721</td>
      <td>BeatAML</td>
      <td>0.7264</td>
      <td>NaN</td>
      <td>7.525</td>
      <td>3.249</td>
      <td>0.0277</td>
      <td>0.5343</td>
      <td>4.000</td>
      <td>0.4609</td>
      <td>0.5391</td>
      <td>0.4010</td>
    </tr>
    <tr>
      <th>4</th>
      <td>synapse</td>
      <td>4634</td>
      <td>SMI_57735</td>
      <td>BeatAML</td>
      <td>0.5588</td>
      <td>NaN</td>
      <td>11.540</td>
      <td>0.000</td>
      <td>-0.0000</td>
      <td>0.5588</td>
      <td>3.611</td>
      <td>0.4412</td>
      <td>0.5588</td>
      <td>0.3791</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>38775</th>
      <td>synapse</td>
      <td>4346</td>
      <td>SMI_414</td>
      <td>BeatAML</td>
      <td>0.4791</td>
      <td>10.94</td>
      <td>11.360</td>
      <td>0.000</td>
      <td>-0.0000</td>
      <td>0.4791</td>
      <td>3.294</td>
      <td>0.5209</td>
      <td>0.4791</td>
      <td>0.4677</td>
    </tr>
    <tr>
      <th>38776</th>
      <td>synapse</td>
      <td>4346</td>
      <td>SMI_1062</td>
      <td>BeatAML</td>
      <td>0.9319</td>
      <td>NaN</td>
      <td>9.805</td>
      <td>0.000</td>
      <td>-0.0000</td>
      <td>0.8638</td>
      <td>0.000</td>
      <td>0.0681</td>
      <td>0.9319</td>
      <td>0.0000</td>
    </tr>
    <tr>
      <th>38777</th>
      <td>synapse</td>
      <td>4346</td>
      <td>SMI_3146</td>
      <td>BeatAML</td>
      <td>0.4195</td>
      <td>10.81</td>
      <td>11.060</td>
      <td>0.000</td>
      <td>-0.0000</td>
      <td>0.4195</td>
      <td>3.247</td>
      <td>0.5805</td>
      <td>0.4195</td>
      <td>0.5339</td>
    </tr>
    <tr>
      <th>38778</th>
      <td>synapse</td>
      <td>4346</td>
      <td>SMI_1137</td>
      <td>BeatAML</td>
      <td>0.8347</td>
      <td>NaN</td>
      <td>0.036</td>
      <td>0.000</td>
      <td>-0.0000</td>
      <td>0.6693</td>
      <td>0.000</td>
      <td>0.1653</td>
      <td>0.8347</td>
      <td>0.0726</td>
    </tr>
    <tr>
      <th>38779</th>
      <td>synapse</td>
      <td>4346</td>
      <td>SMI_181</td>
      <td>BeatAML</td>
      <td>0.7405</td>
      <td>NaN</td>
      <td>11.370</td>
      <td>0.000</td>
      <td>-0.0000</td>
      <td>0.7405</td>
      <td>3.878</td>
      <td>0.2595</td>
      <td>0.7405</td>
      <td>0.1773</td>
    </tr>
  </tbody>
</table>
<p>38780 rows × 14 columns</p>
</div>




```python
merged_df = pd.merge(beataml.experiments[['improve_sample_id', 'improve_drug_id', 'auc']],
                     beataml.drugs[['improve_drug_id', 'isoSMILES']],
                     on='improve_drug_id',
                     how='inner').drop_duplicates()

merged_df = merged_df.dropna(subset=['isoSMILES'])
```


```python
merged_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>improve_sample_id</th>
      <th>improve_drug_id</th>
      <th>auc</th>
      <th>isoSMILES</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>4634</td>
      <td>SMI_43</td>
      <td>0.5716</td>
      <td>C[C@H]1C[C@@H]([C@@H]([C@H](/C=C(/[C@@H]([C@H]...</td>
    </tr>
    <tr>
      <th>284</th>
      <td>4634</td>
      <td>SMI_3222</td>
      <td>0.6310</td>
      <td>CC1=C2C=C(C=CC2=NN1)C3=CC(=CN=C3)OC[C@H](CC4=C...</td>
    </tr>
    <tr>
      <th>569</th>
      <td>4634</td>
      <td>SMI_146</td>
      <td>0.7496</td>
      <td>CN(C)CC[C@H](CSC1=CC=CC=C1)NC2=C(C=C(C=C2)S(=O...</td>
    </tr>
    <tr>
      <th>781</th>
      <td>4634</td>
      <td>SMI_57721</td>
      <td>0.7264</td>
      <td>CC[N+]1=C(N(C2=C1C=C(C=C2)C3=NC4=CC=CC=C4S3)C5...</td>
    </tr>
    <tr>
      <th>846</th>
      <td>4634</td>
      <td>SMI_57735</td>
      <td>0.5588</td>
      <td>C1CCN(CC1)CCOC2=CC=C(C=C2)C3=CN4C(=C(C=N4)C5=C...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>10668667</th>
      <td>4346</td>
      <td>SMI_414</td>
      <td>0.4791</td>
      <td>CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC...</td>
    </tr>
    <tr>
      <th>10668944</th>
      <td>4346</td>
      <td>SMI_1062</td>
      <td>0.9319</td>
      <td>CS(=O)(=O)C1=CC(=C(C=C1)C(=O)NC2=CC(=C(C=C2)Cl...</td>
    </tr>
    <tr>
      <th>10669230</th>
      <td>4346</td>
      <td>SMI_3146</td>
      <td>0.4195</td>
      <td>CC[C@@H]1C(=O)N(C2=CN=C(N=C2N1C(C)C)NC3=C(C=C(...</td>
    </tr>
    <tr>
      <th>10669505</th>
      <td>4346</td>
      <td>SMI_1137</td>
      <td>0.8347</td>
      <td>C1CSCC2=C1N=C(NC2=O)C3=CC=C(C=C3)C(F)(F)F</td>
    </tr>
    <tr>
      <th>10669748</th>
      <td>4346</td>
      <td>SMI_181</td>
      <td>0.7405</td>
      <td>CC1=[N+](C2=C(N1CCOC)C(=O)C3=CC=CC=C3C2=O)CC4=...</td>
    </tr>
  </tbody>
</table>
<p>34850 rows × 4 columns</p>
</div>




```python
# Convert SMILES to Morgan fingerprints
def smiles_to_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)  # Adjust radius and nBits as needed
    fingerprint_array = np.array(fingerprint)
    return fingerprint_array

# smiles_to_fingerprint(merged_df.isoSMILES[0])
# # Apply the function to your SMILES column
merged_df['fingerprint'] = merged_df['isoSMILES'].apply(smiles_to_fingerprint)

fingerprint_map = merged_df[['fingerprint',"improve_drug_id"]]
merged_df = merged_df[["improve_sample_id","improve_drug_id","auc"]]
```


```python
merged_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>improve_sample_id</th>
      <th>improve_drug_id</th>
      <th>auc</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>4634</td>
      <td>SMI_43</td>
      <td>0.5716</td>
    </tr>
    <tr>
      <th>284</th>
      <td>4634</td>
      <td>SMI_3222</td>
      <td>0.6310</td>
    </tr>
    <tr>
      <th>569</th>
      <td>4634</td>
      <td>SMI_146</td>
      <td>0.7496</td>
    </tr>
    <tr>
      <th>781</th>
      <td>4634</td>
      <td>SMI_57721</td>
      <td>0.7264</td>
    </tr>
    <tr>
      <th>846</th>
      <td>4634</td>
      <td>SMI_57735</td>
      <td>0.5588</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>10668667</th>
      <td>4346</td>
      <td>SMI_414</td>
      <td>0.4791</td>
    </tr>
    <tr>
      <th>10668944</th>
      <td>4346</td>
      <td>SMI_1062</td>
      <td>0.9319</td>
    </tr>
    <tr>
      <th>10669230</th>
      <td>4346</td>
      <td>SMI_3146</td>
      <td>0.4195</td>
    </tr>
    <tr>
      <th>10669505</th>
      <td>4346</td>
      <td>SMI_1137</td>
      <td>0.8347</td>
    </tr>
    <tr>
      <th>10669748</th>
      <td>4346</td>
      <td>SMI_181</td>
      <td>0.7405</td>
    </tr>
  </tbody>
</table>
<p>34850 rows × 3 columns</p>
</div>




```python
fingerprint_map
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fingerprint</th>
      <th>improve_drug_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, ...</td>
      <td>SMI_43</td>
    </tr>
    <tr>
      <th>284</th>
      <td>[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
      <td>SMI_3222</td>
    </tr>
    <tr>
      <th>569</th>
      <td>[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, ...</td>
      <td>SMI_146</td>
    </tr>
    <tr>
      <th>781</th>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
      <td>SMI_57721</td>
    </tr>
    <tr>
      <th>846</th>
      <td>[0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, ...</td>
      <td>SMI_57735</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>10668667</th>
      <td>[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, ...</td>
      <td>SMI_414</td>
    </tr>
    <tr>
      <th>10668944</th>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
      <td>SMI_1062</td>
    </tr>
    <tr>
      <th>10669230</th>
      <td>[0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, ...</td>
      <td>SMI_3146</td>
    </tr>
    <tr>
      <th>10669505</th>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
      <td>SMI_1137</td>
    </tr>
    <tr>
      <th>10669748</th>
      <td>[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, ...</td>
      <td>SMI_181</td>
    </tr>
  </tbody>
</table>
<p>34850 rows × 2 columns</p>
</div>




```python
fingerprint_dict = fingerprint_map.set_index('improve_drug_id')['fingerprint'].to_dict()
fingerprint_dict

```




    {'SMI_43': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_3222': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_146': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_57721': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57735': array([0, 0, 1, ..., 0, 0, 0]),
     'SMI_1356': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_949': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_990': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_224': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_274': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_388': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_398': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_98': array([0, 0, 0, ..., 1, 0, 0]),
     'SMI_381': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_145': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_331': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_97': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1302': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_341': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_378': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_225': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_988': array([0, 0, 0, ..., 0, 0, 1]),
     'SMI_57': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_337': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_4071': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1303': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_58': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1003': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_230': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1134': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1397': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1393': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_384': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57737': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_52446': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_3033': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_342': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_95': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57722': array([0, 0, 1, ..., 0, 0, 0]),
     'SMI_52331': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1412': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_1842': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_37': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_57732': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_57724': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57728': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_64': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_338': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_56627': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_41054': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_376': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1060': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_221': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_755': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_125': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57730': array([0, 0, 0, ..., 1, 0, 0]),
     'SMI_3276': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_227': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1057': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_3125': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_318': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_231': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_184': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57723': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1312': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57726': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_148': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_233': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57731': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_53380': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_1108': array([1, 0, 0, ..., 0, 0, 0]),
     'SMI_283': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1403': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1056': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_339': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_63': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_262': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1063': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_276': array([0, 1, 1, ..., 0, 0, 0]),
     'SMI_144': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_55709': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_52447': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_50': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_272': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_232': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57736': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_236': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_223': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_365': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_102': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_51': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_107': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_374': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_369': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_373': array([0, 1, 1, ..., 0, 0, 0]),
     'SMI_188': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_54925': array([0, 0, 1, ..., 0, 0, 0]),
     'SMI_1012': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_995': array([0, 0, 0, ..., 0, 1, 0]),
     'SMI_57733': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_101': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_3148': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1062': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1137': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_181': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_57727': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1963': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_375': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_261': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_278': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1155': array([0, 1, 1, ..., 0, 0, 0]),
     'SMI_1109': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_10': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_52922': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_2582': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_75': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1777': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_1107': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_55184': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_3166': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_410': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_256': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_265': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1327': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_106': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_140': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_1055': array([0, 0, 0, ..., 0, 1, 0]),
     'SMI_60': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_281': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57729': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_1952': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_285': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_414': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_3146': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_1061': array([0, 0, 1, ..., 0, 0, 0]),
     'SMI_1104': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57720': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_349': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_167': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_28': array([1, 1, 0, ..., 0, 0, 0]),
     'SMI_165': array([1, 0, 0, ..., 0, 0, 0]),
     'SMI_1859': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_54398': array([0, 0, 1, ..., 0, 0, 0]),
     'SMI_1984': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_182': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57734': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_53718': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_2064': array([0, 1, 0, ..., 0, 0, 0]),
     'SMI_299': array([0, 0, 0, ..., 0, 0, 0]),
     'SMI_57725': array([0, 0, 0, ..., 0, 0, 0])}




```python

```


```python
# Merge merged_df with transcriptomics_df based on improve_sample_id
tp = pd.merge(beataml.transcriptomics,
                     beataml.proteomics[['improve_sample_id', 'proteomics', 'entrez_id']],
                     on=['improve_sample_id','entrez_id'],
                     how='inner').drop_duplicates()

merged_df = pd.merge(merged_df,
                     tp[['improve_sample_id', 'transcriptomics','proteomics', 'entrez_id']],
                     on='improve_sample_id',
                     how='inner').drop_duplicates()


```


```python
merged_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>improve_sample_id</th>
      <th>improve_drug_id</th>
      <th>auc</th>
      <th>transcriptomics</th>
      <th>proteomics</th>
      <th>entrez_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>7.068854</td>
      <td>-0.135</td>
      <td>8813.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>3.859919</td>
      <td>-0.300</td>
      <td>6359.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>3.859919</td>
      <td>-0.300</td>
      <td>57147.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>9.634914</td>
      <td>1.140</td>
      <td>2268.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>-1.401864</td>
      <td>-0.174</td>
      <td>3075.0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>8148099</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>1.519607</td>
      <td>-0.192</td>
      <td>80755.0</td>
    </tr>
    <tr>
      <th>8148100</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>2.088872</td>
      <td>-0.260</td>
      <td>55957.0</td>
    </tr>
    <tr>
      <th>8148101</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>2.154935</td>
      <td>0.865</td>
      <td>6023.0</td>
    </tr>
    <tr>
      <th>8148102</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>7.165478</td>
      <td>0.472</td>
      <td>11165.0</td>
    </tr>
    <tr>
      <th>8148103</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>-0.435228</td>
      <td>0.365</td>
      <td>6047.0</td>
    </tr>
  </tbody>
</table>
<p>8148104 rows × 6 columns</p>
</div>




```python
#Impute missing proteomics (and transcriptomics if there are any) based on global mean. 
columns_to_fill = ['transcriptomics', 'proteomics']

for i in columns_to_fill:
    if i in merged_df.columns[merged_df.isnull().any(axis=0)]:
        merged_df[i].fillna(merged_df[i].mean(), inplace=True)

```


```python
merged_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>improve_sample_id</th>
      <th>improve_drug_id</th>
      <th>auc</th>
      <th>transcriptomics</th>
      <th>proteomics</th>
      <th>entrez_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>7.068854</td>
      <td>-0.135</td>
      <td>8813.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>3.859919</td>
      <td>-0.300</td>
      <td>6359.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>3.859919</td>
      <td>-0.300</td>
      <td>57147.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>9.634914</td>
      <td>1.140</td>
      <td>2268.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>-1.401864</td>
      <td>-0.174</td>
      <td>3075.0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>8148099</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>1.519607</td>
      <td>-0.192</td>
      <td>80755.0</td>
    </tr>
    <tr>
      <th>8148100</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>2.088872</td>
      <td>-0.260</td>
      <td>55957.0</td>
    </tr>
    <tr>
      <th>8148101</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>2.154935</td>
      <td>0.865</td>
      <td>6023.0</td>
    </tr>
    <tr>
      <th>8148102</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>7.165478</td>
      <td>0.472</td>
      <td>11165.0</td>
    </tr>
    <tr>
      <th>8148103</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>-0.435228</td>
      <td>0.365</td>
      <td>6047.0</td>
    </tr>
  </tbody>
</table>
<p>8148104 rows × 6 columns</p>
</div>




```python
merged_df['fingerprint'] = merged_df['improve_drug_id'].map(fingerprint_dict)
# merged_df = merged_df[["improve_sample_id"
```


```python

```


```python
merged_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>improve_sample_id</th>
      <th>improve_drug_id</th>
      <th>auc</th>
      <th>transcriptomics</th>
      <th>proteomics</th>
      <th>entrez_id</th>
      <th>fingerprint</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>7.068854</td>
      <td>-0.135</td>
      <td>8813.0</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, ...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>3.859919</td>
      <td>-0.300</td>
      <td>6359.0</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, ...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>3.859919</td>
      <td>-0.300</td>
      <td>57147.0</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, ...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>9.634914</td>
      <td>1.140</td>
      <td>2268.0</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, ...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>3986</td>
      <td>SMI_236</td>
      <td>0.0000</td>
      <td>-1.401864</td>
      <td>-0.174</td>
      <td>3075.0</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, ...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>8148099</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>1.519607</td>
      <td>-0.192</td>
      <td>80755.0</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>8148100</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>2.088872</td>
      <td>-0.260</td>
      <td>55957.0</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>8148101</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>2.154935</td>
      <td>0.865</td>
      <td>6023.0</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>8148102</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>7.165478</td>
      <td>0.472</td>
      <td>11165.0</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>8148103</th>
      <td>4007</td>
      <td>SMI_57727</td>
      <td>0.9307</td>
      <td>-0.435228</td>
      <td>0.365</td>
      <td>6047.0</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
  </tbody>
</table>
<p>8148104 rows × 7 columns</p>
</div>




```python
scaler = MinMaxScaler()
merged_df[["transcriptomics","proteomics"]] = scaler.fit_transform(merged_df[["transcriptomics","proteomics"]])

```


```python

```


```python
#Create A Hashmap of {improve_sample_id: {improve_drug_id: { entrez_id: transcriptomics, proteomics, auc]}}

data_dict = {}

# Iterate over the DataFrame rows to populate the nested dictionary
for index, row in merged_df.iterrows():
    improve_sample_id = row['improve_sample_id']
    improve_drug_id = row['improve_drug_id']
    auc = row['auc']
    transcriptomics = row['transcriptomics']
    proteomics = row['proteomics']
    entrez_id = row['entrez_id']
    # print(improve_sample_id,improve_drug_id,auc,transcriptomics,proteomics,entrez_id)

    # Initialize improve_sample_id key if not present
    if improve_sample_id not in data_dict:
        data_dict[improve_sample_id] = {}
    # print(data_dict)
#     # Initialize improve_drug_id key if not present
    if improve_drug_id not in data_dict[improve_sample_id]:
        data_dict[improve_sample_id][improve_drug_id] = {}
#     # Append data to the nested dictionary
    if entrez_id not in data_dict[improve_sample_id][improve_drug_id]:
        data_dict[improve_sample_id][improve_drug_id][entrez_id] = {
            'transcriptomics': int,
            'proteomics': int,
            'auc': auc
        }
    data_dict[improve_sample_id][improve_drug_id][entrez_id]['transcriptomics'] = transcriptomics 
    data_dict[improve_sample_id][improve_drug_id][entrez_id]['proteomics'] = proteomics 



```


```python
#Write [transcriptomics], [proteomics], and auc data to aggregated_df. These are all written in the same order as the unique_entrez_ids set.

# # Extract the list of unique entrez_ids
unique_entrez_ids = set(merged_df.entrez_id.unique())

# # Construct the final DataFrame
final_data = []
for improve_sample_id, sample_id_data in data_dict.items():
    for improve_drug_id, drug_id_data in sample_id_data.items():
        transcriptomics = []
        proteomics = []
        for entrez_id in unique_entrez_ids:
            if entrez_id in drug_id_data:
                transcriptomics.append(drug_id_data[entrez_id]['transcriptomics'])
                proteomics.append(drug_id_data[entrez_id]['proteomics'])
                auc = drug_id_data[entrez_id]['auc']
        final_data.append([auc, transcriptomics, proteomics, improve_sample_id,improve_drug_id])

# # Create the DataFrame
aggregated_df = pd.DataFrame(final_data, columns=['auc', 'transcriptomics', 'proteomics', 'improve_sample_id', "improve_drug_id"])
aggregated_df


```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>auc</th>
      <th>transcriptomics</th>
      <th>proteomics</th>
      <th>improve_sample_id</th>
      <th>improve_drug_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.0000</td>
      <td>[0.23608062994353704, 0.5037978202291047, 0.48...</td>
      <td>[0.3384677419354839, 0.37193548387096775, 0.38...</td>
      <td>3986</td>
      <td>SMI_236</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.7431</td>
      <td>[0.23608062994353704, 0.5037978202291047, 0.48...</td>
      <td>[0.3384677419354839, 0.37193548387096775, 0.38...</td>
      <td>3986</td>
      <td>SMI_223</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0.6135</td>
      <td>[0.23608062994353704, 0.5037978202291047, 0.48...</td>
      <td>[0.3384677419354839, 0.37193548387096775, 0.38...</td>
      <td>3986</td>
      <td>SMI_365</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.7425</td>
      <td>[0.18000067030285077, 0.4920917394967491, 0.47...</td>
      <td>[0.43870967741935485, 0.3573387096774194, 0.38...</td>
      <td>3987</td>
      <td>SMI_54925</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.0000</td>
      <td>[0.18000067030285077, 0.4920917394967491, 0.47...</td>
      <td>[0.43870967741935485, 0.3573387096774194, 0.38...</td>
      <td>3987</td>
      <td>SMI_1012</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>988</th>
      <td>0.8771</td>
      <td>[0.29431256693950825, 0.5401532345082279, 0.47...</td>
      <td>[0.34040322580645166, 0.3787903225806452, 0.37...</td>
      <td>4007</td>
      <td>SMI_3148</td>
    </tr>
    <tr>
      <th>989</th>
      <td>1.0000</td>
      <td>[0.29431256693950825, 0.5401532345082279, 0.47...</td>
      <td>[0.34040322580645166, 0.3787903225806452, 0.37...</td>
      <td>4007</td>
      <td>SMI_1062</td>
    </tr>
    <tr>
      <th>990</th>
      <td>0.5961</td>
      <td>[0.29431256693950825, 0.5401532345082279, 0.47...</td>
      <td>[0.34040322580645166, 0.3787903225806452, 0.37...</td>
      <td>4007</td>
      <td>SMI_1137</td>
    </tr>
    <tr>
      <th>991</th>
      <td>0.9997</td>
      <td>[0.29431256693950825, 0.5401532345082279, 0.47...</td>
      <td>[0.34040322580645166, 0.3787903225806452, 0.37...</td>
      <td>4007</td>
      <td>SMI_181</td>
    </tr>
    <tr>
      <th>992</th>
      <td>0.9307</td>
      <td>[0.29431256693950825, 0.5401532345082279, 0.47...</td>
      <td>[0.34040322580645166, 0.3787903225806452, 0.37...</td>
      <td>4007</td>
      <td>SMI_57727</td>
    </tr>
  </tbody>
</table>
<p>993 rows × 5 columns</p>
</div>




```python
# Map in the fingerprint by improve_drug_id. Remove improve_drug_id.
aggregated_df['fingerprint'] = aggregated_df['improve_drug_id'].map(fingerprint_dict)
aggregated_df = aggregated_df[["auc","transcriptomics","proteomics","fingerprint"]]
```


```python
# Training/Test/Validation Data is all under the name training_data.
training_data = aggregated_df
training_data

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>auc</th>
      <th>transcriptomics</th>
      <th>proteomics</th>
      <th>fingerprint</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.0000</td>
      <td>[0.23608062994353704, 0.5037978202291047, 0.48...</td>
      <td>[0.3384677419354839, 0.37193548387096775, 0.38...</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, ...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.7431</td>
      <td>[0.23608062994353704, 0.5037978202291047, 0.48...</td>
      <td>[0.3384677419354839, 0.37193548387096775, 0.38...</td>
      <td>[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, ...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0.6135</td>
      <td>[0.23608062994353704, 0.5037978202291047, 0.48...</td>
      <td>[0.3384677419354839, 0.37193548387096775, 0.38...</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, ...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.7425</td>
      <td>[0.18000067030285077, 0.4920917394967491, 0.47...</td>
      <td>[0.43870967741935485, 0.3573387096774194, 0.38...</td>
      <td>[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.0000</td>
      <td>[0.18000067030285077, 0.4920917394967491, 0.47...</td>
      <td>[0.43870967741935485, 0.3573387096774194, 0.38...</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, ...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>988</th>
      <td>0.8771</td>
      <td>[0.29431256693950825, 0.5401532345082279, 0.47...</td>
      <td>[0.34040322580645166, 0.3787903225806452, 0.37...</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>989</th>
      <td>1.0000</td>
      <td>[0.29431256693950825, 0.5401532345082279, 0.47...</td>
      <td>[0.34040322580645166, 0.3787903225806452, 0.37...</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>990</th>
      <td>0.5961</td>
      <td>[0.29431256693950825, 0.5401532345082279, 0.47...</td>
      <td>[0.34040322580645166, 0.3787903225806452, 0.37...</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>991</th>
      <td>0.9997</td>
      <td>[0.29431256693950825, 0.5401532345082279, 0.47...</td>
      <td>[0.34040322580645166, 0.3787903225806452, 0.37...</td>
      <td>[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, ...</td>
    </tr>
    <tr>
      <th>992</th>
      <td>0.9307</td>
      <td>[0.29431256693950825, 0.5401532345082279, 0.47...</td>
      <td>[0.34040322580645166, 0.3787903225806452, 0.37...</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
  </tbody>
</table>
<p>993 rows × 4 columns</p>
</div>




```python

# Splitting the merged_data DataFrame into features (X) and target variable (y)
X = training_data[['transcriptomics', 'proteomics',"fingerprint"]]  # These columns are your features
y = training_data['auc']  # Using 'auc' column as the target variable

```


```python

# Splitting data into training, validation, and testing sets
X_train_val, X_test, y_train_val, y_test = train_test_split(X, y, test_size=0.1, random_state=42)
X_train, X_val, y_train, y_val = train_test_split(X_train_val, y_train_val, test_size=0.1111, random_state=42) # 0.1111 * 90 = 10%

```


```python
#These are the array lengths of each feature set.
num_proteomics = 7284
num_transcriptomics = 7284  
num_fingerprint = 1024  

#Define inputs
transcriptomics_input = keras.Input(
    shape=(num_transcriptomics,), name="transcriptomics"
)  

proteomics_input = keras.Input(
    shape=(num_proteomics,), name="proteomics"
)  

fingerprint_input = keras.Input(
    shape=(num_fingerprint,), name="fingerprint"
)  

# You may adjust these layers and activation functions to get better (or worse) results.
# Merge all available features into a single large vector via concatenation
x = layers.concatenate([transcriptomics_input, proteomics_input, fingerprint_input])
x = layers.Dense(16, activation='relu')(x)
x = layers.Dense(64, activation="relu")(x)
x = layers.Dropout(0.5)(x)
x = layers.Dense(32, activation='relu')(x)
x = layers.Dense(12, activation='relu')(x)

# Priority prediction is here. You may add layers (and an output) after this too.
# This just allows for an early prediction in case you have a large datset and want preliminary results.
priority_pred = layers.Dense(1, name="priority",activation='relu')(x)

# Instantiate an end-to-end model predicting both priority and department
model = keras.Model(
    inputs=[transcriptomics_input, proteomics_input, fingerprint_input],
    outputs={"priority": priority_pred},
)
```


```python
keras.utils.plot_model(model, "multi_input_model.png", show_shapes=True)

```

    You must install pydot (`pip install pydot`) and install graphviz (see instructions at https://graphviz.gitlab.io/download/) for plot_model to work.



```python
model.compile(
    optimizer=keras.optimizers.Adam(learning_rate=0.0001),
    loss={
        "priority": keras.losses.MeanSquaredError(),
    },
    metrics=[MeanAbsoluteError()],
    loss_weights={"priority": 1},
)

```


```python
history = model.fit(
    [np.array(X_train['transcriptomics'].tolist()),
     np.array(X_train['proteomics'].tolist()),
     np.array(X_train['fingerprint'].tolist())
    ], 
    y_train, 
    epochs=100, 
    batch_size=32, 
    validation_data=(
        [
            np.array(X_val['transcriptomics'].tolist()), 
            np.array(X_val['proteomics'].tolist()), 
            np.array(X_val['fingerprint'].tolist())
        ], 
        y_val)
)

```

    Epoch 1/100


    2024-03-07 22:06:39.578569: W tensorflow/tsl/platform/profile_utils/cpu_utils.cc:128] Failed to get CPU frequency: 0 Hz


    25/25 [==============================] - 1s 16ms/step - loss: 0.3187 - mean_absolute_error: 0.4867 - val_loss: 0.2245 - val_mean_absolute_error: 0.4414
    Epoch 2/100
    25/25 [==============================] - 0s 5ms/step - loss: 0.2099 - mean_absolute_error: 0.3757 - val_loss: 0.1138 - val_mean_absolute_error: 0.3021
    Epoch 3/100
    25/25 [==============================] - 0s 5ms/step - loss: 0.1783 - mean_absolute_error: 0.3399 - val_loss: 0.0712 - val_mean_absolute_error: 0.2290
    Epoch 4/100
    25/25 [==============================] - 0s 4ms/step - loss: 0.1419 - mean_absolute_error: 0.3088 - val_loss: 0.0421 - val_mean_absolute_error: 0.1584
    Epoch 5/100
    25/25 [==============================] - 0s 5ms/step - loss: 0.1281 - mean_absolute_error: 0.2879 - val_loss: 0.0451 - val_mean_absolute_error: 0.1690
    Epoch 6/100
    25/25 [==============================] - 0s 5ms/step - loss: 0.1105 - mean_absolute_error: 0.2707 - val_loss: 0.0426 - val_mean_absolute_error: 0.1616
    Epoch 7/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0958 - mean_absolute_error: 0.2483 - val_loss: 0.0388 - val_mean_absolute_error: 0.1455
    Epoch 8/100
    25/25 [==============================] - 0s 4ms/step - loss: 0.0988 - mean_absolute_error: 0.2488 - val_loss: 0.0382 - val_mean_absolute_error: 0.1433
    Epoch 9/100
    25/25 [==============================] - 0s 4ms/step - loss: 0.0867 - mean_absolute_error: 0.2319 - val_loss: 0.0379 - val_mean_absolute_error: 0.1409
    Epoch 10/100
    25/25 [==============================] - 0s 5ms/step - loss: 0.0809 - mean_absolute_error: 0.2240 - val_loss: 0.0400 - val_mean_absolute_error: 0.1426
    Epoch 11/100
    25/25 [==============================] - 0s 4ms/step - loss: 0.0820 - mean_absolute_error: 0.2302 - val_loss: 0.0414 - val_mean_absolute_error: 0.1450
    Epoch 12/100
    25/25 [==============================] - 0s 7ms/step - loss: 0.0767 - mean_absolute_error: 0.2233 - val_loss: 0.0391 - val_mean_absolute_error: 0.1558
    Epoch 13/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0804 - mean_absolute_error: 0.2246 - val_loss: 0.0366 - val_mean_absolute_error: 0.1444
    Epoch 14/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0736 - mean_absolute_error: 0.2140 - val_loss: 0.0364 - val_mean_absolute_error: 0.1438
    Epoch 15/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0742 - mean_absolute_error: 0.2169 - val_loss: 0.0420 - val_mean_absolute_error: 0.1664
    Epoch 16/100
    25/25 [==============================] - 0s 4ms/step - loss: 0.0803 - mean_absolute_error: 0.2277 - val_loss: 0.0417 - val_mean_absolute_error: 0.1655
    Epoch 17/100
    25/25 [==============================] - 0s 6ms/step - loss: 0.0665 - mean_absolute_error: 0.2101 - val_loss: 0.0380 - val_mean_absolute_error: 0.1364
    Epoch 18/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0641 - mean_absolute_error: 0.2001 - val_loss: 0.0376 - val_mean_absolute_error: 0.1515
    Epoch 19/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0650 - mean_absolute_error: 0.2043 - val_loss: 0.0381 - val_mean_absolute_error: 0.1360
    Epoch 20/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0579 - mean_absolute_error: 0.1924 - val_loss: 0.0430 - val_mean_absolute_error: 0.1476
    Epoch 21/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0610 - mean_absolute_error: 0.1939 - val_loss: 0.0606 - val_mean_absolute_error: 0.2099
    Epoch 22/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0576 - mean_absolute_error: 0.1909 - val_loss: 0.0352 - val_mean_absolute_error: 0.1401
    Epoch 23/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0554 - mean_absolute_error: 0.1864 - val_loss: 0.0348 - val_mean_absolute_error: 0.1357
    Epoch 24/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0578 - mean_absolute_error: 0.1866 - val_loss: 0.0401 - val_mean_absolute_error: 0.1606
    Epoch 25/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0553 - mean_absolute_error: 0.1847 - val_loss: 0.0380 - val_mean_absolute_error: 0.1538
    Epoch 26/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0481 - mean_absolute_error: 0.1711 - val_loss: 0.0354 - val_mean_absolute_error: 0.1427
    Epoch 27/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0528 - mean_absolute_error: 0.1798 - val_loss: 0.0383 - val_mean_absolute_error: 0.1354
    Epoch 28/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0516 - mean_absolute_error: 0.1753 - val_loss: 0.0355 - val_mean_absolute_error: 0.1435
    Epoch 29/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0506 - mean_absolute_error: 0.1744 - val_loss: 0.0349 - val_mean_absolute_error: 0.1393
    Epoch 30/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0567 - mean_absolute_error: 0.1873 - val_loss: 0.0591 - val_mean_absolute_error: 0.2049
    Epoch 31/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0544 - mean_absolute_error: 0.1824 - val_loss: 0.0346 - val_mean_absolute_error: 0.1352
    Epoch 32/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0459 - mean_absolute_error: 0.1687 - val_loss: 0.0368 - val_mean_absolute_error: 0.1491
    Epoch 33/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0447 - mean_absolute_error: 0.1643 - val_loss: 0.0348 - val_mean_absolute_error: 0.1342
    Epoch 34/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0425 - mean_absolute_error: 0.1626 - val_loss: 0.0347 - val_mean_absolute_error: 0.1354
    Epoch 35/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0422 - mean_absolute_error: 0.1596 - val_loss: 0.0381 - val_mean_absolute_error: 0.1530
    Epoch 36/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0454 - mean_absolute_error: 0.1644 - val_loss: 0.0350 - val_mean_absolute_error: 0.1382
    Epoch 37/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0432 - mean_absolute_error: 0.1594 - val_loss: 0.0375 - val_mean_absolute_error: 0.1358
    Epoch 38/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0435 - mean_absolute_error: 0.1617 - val_loss: 0.0348 - val_mean_absolute_error: 0.1345
    Epoch 39/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0390 - mean_absolute_error: 0.1518 - val_loss: 0.0349 - val_mean_absolute_error: 0.1369
    Epoch 40/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0416 - mean_absolute_error: 0.1586 - val_loss: 0.0365 - val_mean_absolute_error: 0.1353
    Epoch 41/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0380 - mean_absolute_error: 0.1496 - val_loss: 0.0352 - val_mean_absolute_error: 0.1354
    Epoch 42/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0406 - mean_absolute_error: 0.1561 - val_loss: 0.0373 - val_mean_absolute_error: 0.1367
    Epoch 43/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0365 - mean_absolute_error: 0.1461 - val_loss: 0.0351 - val_mean_absolute_error: 0.1364
    Epoch 44/100
    25/25 [==============================] - 0s 6ms/step - loss: 0.0397 - mean_absolute_error: 0.1514 - val_loss: 0.0364 - val_mean_absolute_error: 0.1448
    Epoch 45/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0391 - mean_absolute_error: 0.1487 - val_loss: 0.0363 - val_mean_absolute_error: 0.1364
    Epoch 46/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0402 - mean_absolute_error: 0.1519 - val_loss: 0.0361 - val_mean_absolute_error: 0.1366
    Epoch 47/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0392 - mean_absolute_error: 0.1544 - val_loss: 0.0362 - val_mean_absolute_error: 0.1411
    Epoch 48/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0412 - mean_absolute_error: 0.1554 - val_loss: 0.0357 - val_mean_absolute_error: 0.1373
    Epoch 49/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0386 - mean_absolute_error: 0.1504 - val_loss: 0.0372 - val_mean_absolute_error: 0.1377
    Epoch 50/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0357 - mean_absolute_error: 0.1443 - val_loss: 0.0370 - val_mean_absolute_error: 0.1373
    Epoch 51/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0387 - mean_absolute_error: 0.1513 - val_loss: 0.0429 - val_mean_absolute_error: 0.1456
    Epoch 52/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0415 - mean_absolute_error: 0.1567 - val_loss: 0.0508 - val_mean_absolute_error: 0.1599
    Epoch 53/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0410 - mean_absolute_error: 0.1550 - val_loss: 0.0362 - val_mean_absolute_error: 0.1399
    Epoch 54/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0349 - mean_absolute_error: 0.1401 - val_loss: 0.0371 - val_mean_absolute_error: 0.1441
    Epoch 55/100
    25/25 [==============================] - 0s 4ms/step - loss: 0.0441 - mean_absolute_error: 0.1634 - val_loss: 0.0546 - val_mean_absolute_error: 0.1954
    Epoch 56/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0404 - mean_absolute_error: 0.1574 - val_loss: 0.0385 - val_mean_absolute_error: 0.1388
    Epoch 57/100
    25/25 [==============================] - 0s 7ms/step - loss: 0.0332 - mean_absolute_error: 0.1395 - val_loss: 0.0406 - val_mean_absolute_error: 0.1418
    Epoch 58/100
    25/25 [==============================] - 0s 4ms/step - loss: 0.0370 - mean_absolute_error: 0.1474 - val_loss: 0.0523 - val_mean_absolute_error: 0.1635
    Epoch 59/100
    25/25 [==============================] - 0s 4ms/step - loss: 0.0356 - mean_absolute_error: 0.1432 - val_loss: 0.0447 - val_mean_absolute_error: 0.1492
    Epoch 60/100
    25/25 [==============================] - 0s 5ms/step - loss: 0.0360 - mean_absolute_error: 0.1446 - val_loss: 0.0370 - val_mean_absolute_error: 0.1425
    Epoch 61/100
    25/25 [==============================] - 0s 4ms/step - loss: 0.0357 - mean_absolute_error: 0.1438 - val_loss: 0.0379 - val_mean_absolute_error: 0.1387
    Epoch 62/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0361 - mean_absolute_error: 0.1469 - val_loss: 0.0368 - val_mean_absolute_error: 0.1416
    Epoch 63/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0341 - mean_absolute_error: 0.1407 - val_loss: 0.0395 - val_mean_absolute_error: 0.1411
    Epoch 64/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0364 - mean_absolute_error: 0.1440 - val_loss: 0.0432 - val_mean_absolute_error: 0.1468
    Epoch 65/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0333 - mean_absolute_error: 0.1380 - val_loss: 0.0367 - val_mean_absolute_error: 0.1389
    Epoch 66/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0346 - mean_absolute_error: 0.1425 - val_loss: 0.0389 - val_mean_absolute_error: 0.1404
    Epoch 67/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0374 - mean_absolute_error: 0.1482 - val_loss: 0.0396 - val_mean_absolute_error: 0.1413
    Epoch 68/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0340 - mean_absolute_error: 0.1386 - val_loss: 0.0366 - val_mean_absolute_error: 0.1403
    Epoch 69/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0351 - mean_absolute_error: 0.1398 - val_loss: 0.0378 - val_mean_absolute_error: 0.1390
    Epoch 70/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0330 - mean_absolute_error: 0.1348 - val_loss: 0.0376 - val_mean_absolute_error: 0.1392
    Epoch 71/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0339 - mean_absolute_error: 0.1389 - val_loss: 0.0383 - val_mean_absolute_error: 0.1402
    Epoch 72/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0322 - mean_absolute_error: 0.1342 - val_loss: 0.0375 - val_mean_absolute_error: 0.1396
    Epoch 73/100
    25/25 [==============================] - 0s 5ms/step - loss: 0.0333 - mean_absolute_error: 0.1385 - val_loss: 0.0410 - val_mean_absolute_error: 0.1439
    Epoch 74/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0312 - mean_absolute_error: 0.1330 - val_loss: 0.0374 - val_mean_absolute_error: 0.1432
    Epoch 75/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0328 - mean_absolute_error: 0.1367 - val_loss: 0.0448 - val_mean_absolute_error: 0.1494
    Epoch 76/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0318 - mean_absolute_error: 0.1333 - val_loss: 0.0399 - val_mean_absolute_error: 0.1425
    Epoch 77/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0304 - mean_absolute_error: 0.1323 - val_loss: 0.0549 - val_mean_absolute_error: 0.1687
    Epoch 78/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0343 - mean_absolute_error: 0.1414 - val_loss: 0.0498 - val_mean_absolute_error: 0.1582
    Epoch 79/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0340 - mean_absolute_error: 0.1409 - val_loss: 0.0525 - val_mean_absolute_error: 0.1639
    Epoch 80/100
    25/25 [==============================] - 0s 4ms/step - loss: 0.0377 - mean_absolute_error: 0.1504 - val_loss: 0.0481 - val_mean_absolute_error: 0.1555
    Epoch 81/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0335 - mean_absolute_error: 0.1385 - val_loss: 0.0416 - val_mean_absolute_error: 0.1448
    Epoch 82/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0305 - mean_absolute_error: 0.1309 - val_loss: 0.0388 - val_mean_absolute_error: 0.1411
    Epoch 83/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0328 - mean_absolute_error: 0.1383 - val_loss: 0.0377 - val_mean_absolute_error: 0.1402
    Epoch 84/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0319 - mean_absolute_error: 0.1346 - val_loss: 0.0377 - val_mean_absolute_error: 0.1406
    Epoch 85/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0313 - mean_absolute_error: 0.1343 - val_loss: 0.0402 - val_mean_absolute_error: 0.1432
    Epoch 86/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0295 - mean_absolute_error: 0.1294 - val_loss: 0.0370 - val_mean_absolute_error: 0.1422
    Epoch 87/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0299 - mean_absolute_error: 0.1296 - val_loss: 0.0372 - val_mean_absolute_error: 0.1411
    Epoch 88/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0315 - mean_absolute_error: 0.1334 - val_loss: 0.0368 - val_mean_absolute_error: 0.1437
    Epoch 89/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0320 - mean_absolute_error: 0.1355 - val_loss: 0.0401 - val_mean_absolute_error: 0.1433
    Epoch 90/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0295 - mean_absolute_error: 0.1295 - val_loss: 0.0377 - val_mean_absolute_error: 0.1409
    Epoch 91/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0287 - mean_absolute_error: 0.1268 - val_loss: 0.0381 - val_mean_absolute_error: 0.1412
    Epoch 92/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0321 - mean_absolute_error: 0.1355 - val_loss: 0.0456 - val_mean_absolute_error: 0.1517
    Epoch 93/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0311 - mean_absolute_error: 0.1331 - val_loss: 0.0415 - val_mean_absolute_error: 0.1456
    Epoch 94/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0310 - mean_absolute_error: 0.1324 - val_loss: 0.0411 - val_mean_absolute_error: 0.1450
    Epoch 95/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0296 - mean_absolute_error: 0.1299 - val_loss: 0.0429 - val_mean_absolute_error: 0.1474
    Epoch 96/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0299 - mean_absolute_error: 0.1307 - val_loss: 0.0374 - val_mean_absolute_error: 0.1414
    Epoch 97/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0296 - mean_absolute_error: 0.1288 - val_loss: 0.0449 - val_mean_absolute_error: 0.1504
    Epoch 98/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0317 - mean_absolute_error: 0.1353 - val_loss: 0.0381 - val_mean_absolute_error: 0.1505
    Epoch 99/100
    25/25 [==============================] - 0s 3ms/step - loss: 0.0303 - mean_absolute_error: 0.1298 - val_loss: 0.0413 - val_mean_absolute_error: 0.1455
    Epoch 100/100
    25/25 [==============================] - 0s 2ms/step - loss: 0.0305 - mean_absolute_error: 0.1282 - val_loss: 0.0408 - val_mean_absolute_error: 0.1447



```python
losses = model.evaluate(
    [np.array(X_test['transcriptomics'].tolist()),
     np.array(X_test['proteomics'].tolist()),
     np.array(X_test['fingerprint'].tolist())
    ], 
    y_test, 
    return_dict=True) 
print(losses)
```

    4/4 [==============================] - 0s 1ms/step - loss: 0.0218 - mean_absolute_error: 0.1116
    {'loss': 0.021844929084181786, 'mean_absolute_error': 0.11160688102245331}



```python

```


```python
# summarize history for accuracy
plt.plot(history.history['val_mean_absolute_error'])
plt.plot(history.history['mean_absolute_error'])
plt.title('Mean Absolute Error of Deep Learning Model')
plt.ylabel('Mean Absolute Error')
plt.xlabel('epoch')
plt.legend(['Validation Set','Training Set',], loc='upper right')
plt.show()
# summarize history for loss
plt.plot(history.history['val_loss'])
plt.plot(history.history['loss'])
plt.title('Loss Function of Deep Learning Model')
plt.ylabel('Loss')
plt.xlabel('epoch')
plt.legend(['Validation Set','Training Set'], loc='upper right')
plt.show()
```


    
![png](assets/markdown/output_34_0.png)



    
![png](assets/markdown/output_34_1.png)
    



```python

# make a prediction
predict = model.predict([np.array(X_test['transcriptomics'].tolist()),
     np.array(X_test['proteomics'].tolist()),
     np.array(X_test['fingerprint'].tolist())
    ])
```

    4/4 [==============================] - 0s 1ms/step



```python
new_df = pd.DataFrame({
    'Predicted Values': predict["priority"].tolist(),
    'True Values': y_test
})
new_df['Predicted Values'] = new_df['Predicted Values'].apply(lambda x: x[0] if isinstance(x, list) and len(x) == 1 else x)
sorted_df = new_df.sort_values(by='True Values', ascending=True)  # Change 'ascending' to False for descending order
sorted_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Predicted Values</th>
      <th>True Values</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>830</th>
      <td>0.689598</td>
      <td>0.0000</td>
    </tr>
    <tr>
      <th>198</th>
      <td>0.613777</td>
      <td>0.3066</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.576904</td>
      <td>0.3380</td>
    </tr>
    <tr>
      <th>355</th>
      <td>0.569279</td>
      <td>0.3572</td>
    </tr>
    <tr>
      <th>656</th>
      <td>0.576972</td>
      <td>0.4010</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>485</th>
      <td>0.907541</td>
      <td>0.9939</td>
    </tr>
    <tr>
      <th>371</th>
      <td>0.961458</td>
      <td>0.9997</td>
    </tr>
    <tr>
      <th>323</th>
      <td>0.957583</td>
      <td>1.0000</td>
    </tr>
    <tr>
      <th>430</th>
      <td>0.911885</td>
      <td>1.0000</td>
    </tr>
    <tr>
      <th>362</th>
      <td>0.750975</td>
      <td>1.0000</td>
    </tr>
  </tbody>
</table>
<p>100 rows × 2 columns</p>
</div>




```python
true_values = np.array(sorted_df["True Values"])
predicted_values = np.array(sorted_df["Predicted Values"])

# Plot the true values and predicted values
plt.plot(true_values)
plt.plot(predicted_values)

# Fit trendlines
true_values_trendline = LinearRegression().fit(np.arange(len(true_values)).reshape(-1, 1), true_values.reshape(-1, 1)).predict(np.arange(len(true_values)).reshape(-1, 1))
predicted_values_trendline = LinearRegression().fit(np.arange(len(predicted_values)).reshape(-1, 1), predicted_values.reshape(-1, 1)).predict(np.arange(len(predicted_values)).reshape(-1, 1))

# Plot trendlines
plt.plot(true_values_trendline, linestyle='--', color='blue', alpha=0.3)
plt.plot(predicted_values_trendline, linestyle='--', color='orange', alpha=0.3)

# Customize plot
plt.title('Predicted AUC vs True AUC')
plt.ylabel('Drug Response AUC')
plt.xlabel('Values')
plt.legend(['True Values', 'Predicted Values', 'True Values Trendline', 'Predicted Values Trendline'], loc='lower right')
plt.show()
```


    
![png](assets/markdown/output_37_0.png)
    



```python
true_values = np.array(new_df["True Values"])
predicted_values = np.array(new_df["Predicted Values"])

# Plot the true values and predicted values
plt.plot(true_values)
plt.plot(predicted_values)

# Fit trendlines
true_values_trendline = LinearRegression().fit(np.arange(len(true_values)).reshape(-1, 1), true_values.reshape(-1, 1)).predict(np.arange(len(true_values)).reshape(-1, 1))
predicted_values_trendline = LinearRegression().fit(np.arange(len(predicted_values)).reshape(-1, 1), predicted_values.reshape(-1, 1)).predict(np.arange(len(predicted_values)).reshape(-1, 1))


# Customize plot
plt.title('Predicted AUC vs True AUC')
plt.ylabel('Drug Response AUC')
plt.xlabel('Values')
plt.legend(['True Values', 'Predicted Values', 'True Values Trendline', 'Predicted Values Trendline'], loc='lower right')
plt.show()
```


    
![png](assets/markdown/output_38_0.png)
    



```python
# Calculate Summary Statistics


# Calculate Mean Absolute Error (MAE)
mae = mean_absolute_error(new_df['True Values'], new_df['Predicted Values'])

# Calculate Root Mean Squared Error (RMSE)
rmse = np.sqrt(mean_squared_error(new_df['True Values'], new_df['Predicted Values']))

# Calculate R-squared (R2) score
r2 = r2_score(new_df['True Values'], new_df['Predicted Values'])

summary_statistics = new_df.describe()

# Print the statistics
print("Mean Absolute Error (MAE):", mae)
print("Root Mean Squared Error (RMSE):", rmse)
print("R-squared (R2) Score:", r2)
print("\n")
# Print summary statistics
print(summary_statistics)
```

    Mean Absolute Error (MAE): 0.11160688844966887
    Root Mean Squared Error (RMSE): 0.14780030079805684
    R-squared (R2) Score: 0.37366931237133216
    
    
           Predicted Values  True Values
    count        100.000000   100.000000
    mean           0.731800     0.723970
    std            0.128604     0.187696
    min            0.477577     0.000000
    25%            0.639343     0.616050
    50%            0.733571     0.736050
    75%            0.838911     0.869575
    max            0.961458     1.000000



```python
#Side by side comparison for first 50 values.
new_df[0:50]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Predicted Values</th>
      <th>True Values</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>918</th>
      <td>0.777755</td>
      <td>0.7367</td>
    </tr>
    <tr>
      <th>525</th>
      <td>0.716599</td>
      <td>0.6711</td>
    </tr>
    <tr>
      <th>567</th>
      <td>0.684600</td>
      <td>0.8426</td>
    </tr>
    <tr>
      <th>656</th>
      <td>0.576972</td>
      <td>0.4010</td>
    </tr>
    <tr>
      <th>915</th>
      <td>0.833463</td>
      <td>0.8296</td>
    </tr>
    <tr>
      <th>429</th>
      <td>0.644341</td>
      <td>0.6900</td>
    </tr>
    <tr>
      <th>855</th>
      <td>0.477577</td>
      <td>0.6126</td>
    </tr>
    <tr>
      <th>711</th>
      <td>0.530543</td>
      <td>0.4982</td>
    </tr>
    <tr>
      <th>174</th>
      <td>0.477577</td>
      <td>0.4917</td>
    </tr>
    <tr>
      <th>604</th>
      <td>0.572155</td>
      <td>0.6209</td>
    </tr>
    <tr>
      <th>865</th>
      <td>0.881899</td>
      <td>0.8678</td>
    </tr>
    <tr>
      <th>449</th>
      <td>0.695015</td>
      <td>0.6786</td>
    </tr>
    <tr>
      <th>777</th>
      <td>0.683768</td>
      <td>0.4146</td>
    </tr>
    <tr>
      <th>580</th>
      <td>0.851610</td>
      <td>0.9087</td>
    </tr>
    <tr>
      <th>76</th>
      <td>0.796451</td>
      <td>0.8499</td>
    </tr>
    <tr>
      <th>371</th>
      <td>0.961458</td>
      <td>0.9997</td>
    </tr>
    <tr>
      <th>884</th>
      <td>0.799178</td>
      <td>0.9071</td>
    </tr>
    <tr>
      <th>136</th>
      <td>0.515090</td>
      <td>0.5763</td>
    </tr>
    <tr>
      <th>158</th>
      <td>0.729075</td>
      <td>0.7231</td>
    </tr>
    <tr>
      <th>290</th>
      <td>0.787732</td>
      <td>0.6100</td>
    </tr>
    <tr>
      <th>673</th>
      <td>0.748977</td>
      <td>0.6931</td>
    </tr>
    <tr>
      <th>321</th>
      <td>0.702108</td>
      <td>0.5028</td>
    </tr>
    <tr>
      <th>757</th>
      <td>0.899316</td>
      <td>0.8691</td>
    </tr>
    <tr>
      <th>70</th>
      <td>0.666315</td>
      <td>0.7865</td>
    </tr>
    <tr>
      <th>355</th>
      <td>0.569279</td>
      <td>0.3572</td>
    </tr>
    <tr>
      <th>359</th>
      <td>0.961239</td>
      <td>0.8362</td>
    </tr>
    <tr>
      <th>107</th>
      <td>0.886554</td>
      <td>0.9438</td>
    </tr>
    <tr>
      <th>265</th>
      <td>0.640555</td>
      <td>0.5975</td>
    </tr>
    <tr>
      <th>825</th>
      <td>0.701618</td>
      <td>0.6212</td>
    </tr>
    <tr>
      <th>139</th>
      <td>0.563990</td>
      <td>0.6810</td>
    </tr>
    <tr>
      <th>184</th>
      <td>0.626397</td>
      <td>0.7691</td>
    </tr>
    <tr>
      <th>708</th>
      <td>0.688292</td>
      <td>0.6542</td>
    </tr>
    <tr>
      <th>622</th>
      <td>0.725480</td>
      <td>0.5147</td>
    </tr>
    <tr>
      <th>941</th>
      <td>0.828256</td>
      <td>0.9047</td>
    </tr>
    <tr>
      <th>305</th>
      <td>0.676587</td>
      <td>0.6286</td>
    </tr>
    <tr>
      <th>809</th>
      <td>0.726858</td>
      <td>0.6427</td>
    </tr>
    <tr>
      <th>306</th>
      <td>0.818704</td>
      <td>0.8063</td>
    </tr>
    <tr>
      <th>767</th>
      <td>0.721371</td>
      <td>0.7445</td>
    </tr>
    <tr>
      <th>863</th>
      <td>0.902120</td>
      <td>0.8096</td>
    </tr>
    <tr>
      <th>323</th>
      <td>0.957583</td>
      <td>1.0000</td>
    </tr>
    <tr>
      <th>59</th>
      <td>0.477577</td>
      <td>0.9167</td>
    </tr>
    <tr>
      <th>298</th>
      <td>0.494210</td>
      <td>0.4433</td>
    </tr>
    <tr>
      <th>668</th>
      <td>0.631534</td>
      <td>0.7201</td>
    </tr>
    <tr>
      <th>617</th>
      <td>0.635709</td>
      <td>0.8713</td>
    </tr>
    <tr>
      <th>370</th>
      <td>0.789941</td>
      <td>0.7637</td>
    </tr>
    <tr>
      <th>23</th>
      <td>0.700036</td>
      <td>0.7924</td>
    </tr>
    <tr>
      <th>30</th>
      <td>0.774137</td>
      <td>0.9460</td>
    </tr>
    <tr>
      <th>816</th>
      <td>0.760774</td>
      <td>0.8721</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.576904</td>
      <td>0.3380</td>
    </tr>
    <tr>
      <th>514</th>
      <td>0.685097</td>
      <td>0.6305</td>
    </tr>
  </tbody>
</table>
</div>


