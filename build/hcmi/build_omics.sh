python 02-getHCMIData.py -m full_manifest.txt -t transcriptomics -o /tmp/hcmi_transcriptomics.csv.gz -g $1 -s $2
python 02-getHCMIData.py -m full_manifest.txt -t copy_number -o /tmp/hcmi_copy_number.csv.gz -g $1 -s $2
python 02-getHCMIData.py -m full_manifest.txt -t mutations -o /tmp/hcmi_mutations.csv.gz -g $1 -s $2

