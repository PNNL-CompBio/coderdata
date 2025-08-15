
#python3 04-experiments-novartispdx.py --token $SYNAPSE_AUTH_TOKEN 

python3 04-experiments-novartis.py --token $SYNAPSE_AUTH_TOKEN -o ~/Projects/CoderData/dev-environment/novartis/novartis_curvedata.tsv
python3 calc_pdx_metrics.py /tmp/novartis_curvedata.tsv --drugfile=/tmp/novartis_drugs.tsv --outprefix=/tmp/novartis --study='Novartis PDX Gao etal 2015' --source='Synapse'