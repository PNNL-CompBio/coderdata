
#python3 04-experiments-novartispdx.py --token $SYNAPSE_AUTH_TOKEN 

python3 -m novartispdx.04-experiments-novartispdx --token $SYNAPSE_AUTH_TOKEN -o ~/Projects/CoderData/dev-environment/novartispdx/novartispdx_curvedata.tsv
python3 utils/calc_pdx_metrics.py /tmp/novartispdx_curvedata.tsv --drugfile=/tmp/novartispdx_drugs.tsv --outprefix=/tmp/novartispdx --study='Novartis PDX Gao etal 2015' --source='Synapse'