Rscript 03_get_drug_response_data.R $SYNAPSE_AUTH_TOKEN $1 $2
/opt/venv/bin/python3 compute_metrics.py /tmp/file.tsv --drugfile=/tmp/mpnstpdx_drugs.tsv
