## HCMI Data

Here we will store the scripts required to process the data from the
[Human Cancer Models
Initiative](https://ocg.cancer.gov/programs/HCMI). 

Currenlty all data collected is part of the [HCMI-CMDC Project on the
GDC](https://portal.gdc.cancer.gov/analysis_page?app=Projects). To
update:

1. Navigate to the [GDC Data
   Portal](https://portal.gdc.cancer.gov/analysis_page?app=Projects),
   and select 'HCMI-CMDC'
2. Click on the 'Cases' button, and select the download button where
   it lists the number of files.
3. This will download the ENTIRE Manifes
4. Filter the manifest for RNASeq, WGS mutations, and copy number
   calls using the following command:
```
 cat ~gdc_manifest.2025-07-08.091940.txt | grep 'mask\|copy\|rna_seq\|md5'
   | grep 'txt\|maf\|tsv\|md5' > new_manifest.txt
 cp new_manifest.txt full_manifest.txt
 
```


Currently the tool require two scripts to build the data:
```
python 01-createHCMISamplesFile.py

python 02-getHCMIData.py -m full_manifest.txt  -t transcriptomics -o transcriptomics.csv

python 02-getHCMIData.py -m full_manifest.txt -t mutations -o mutations.csv

python 02-getHCMIData.py -m full_manifest.txt -t copy_number -o copy_number.csv


```
