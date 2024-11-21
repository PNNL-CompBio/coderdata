Rscript 03_get_drug_response_data.R $SYNAPSE_AUTH_TOKEN $1 $2
/opt/venv/bin/python3 calc_pdx_metrics.py /tmp/curve_data.tsv --drugfile=/tmp/mpnstpdx_drugs.tsv --outprefix=/tmp/mpnstpdx
