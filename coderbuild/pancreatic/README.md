## Pancreatic PDO Data


Here we will store the scripts required to process the omics data from the
Genomic Data Commons together with the drug response data. 

The GDC hosts the panc pdo omcis data, so to update we need an
up-to-date manifest, obtained as follows:


1. Navigate to the [GDC Data
   Portal](https://portal.gdc.cancer.gov/analysis_page?app=Projects),
   and select 'ORGANOID-PANCREATIC'
2. Click on the 'Cases' button, and select the download button where
   it lists the number of files.
3. This will download the ENTIRE Manifest
4. Filter the manifest for RNASeq, WGS mutations, and copy number
   (though i dont think thi dataset has copy number)
   calls using the following command:
```
 cat ~gdc_manifest.2025-07-08.091940.txt | grep 'rna_seq\|md5'
   | 'grep counts\|md5' | grep 'txt\|maf\|tsv\|md5' > new_manifest.txt
 cp new_manifest.txt full_manifest.txt
 
```

The other data is stored [on synapse](https://www.synapse.org/Synapse:syn64597875). 



## Build Docker

```
docker build -f coderbuild/docker/Dockerfile.pancreatic -t pancreatic . 
```

## Run build command
```
python build_dataset.py --dataset pancreatic --build
```

