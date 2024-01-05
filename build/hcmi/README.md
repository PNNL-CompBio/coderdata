## HCMI Data

Here we will store the scripts required to process the data from the [Human Cancer Models Initiative](https://ocg.cancer.gov/programs/HCMI)


Currently the tool require two steps to build the data:
```
python 01-createHCMISamplesFile.py

python 02-getHCMIData.py -m transcriptomics_gdc_manifest.txt  -t transcriptomics -o transcriptomics.csv

python 02-getHCMIData.py -m mutations_manifest_gdc.txt -t mutations -o mutations.csv

python 02-getHCMIData.py -m _manifest.txt -t copy_number -o copy_number.csv


```