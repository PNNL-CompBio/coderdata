## Building cell line data
The cell line data is the first to be built, and requires the
following commands. All scripts write files in to the `/tmp/`
directory, so mounting to that directly will help output the files

```
docker build -f ../dockerfile.cell_line -t cell-line ../
docker run -v $PWD:/tmp/ cell-line Rscript 00-buildGeneFile.R
docker run -v $PWD:/tmp/ cell-line Rscript 01-cellLineSamples.R
docker run -v $PWD:/tmp/ cell-line Rscript 02a-cellLineSanger.R /tmp/genes.csv /tmp/cell_line_samples.csv
docker run -v $PWD:/tmp/ cell-line Rscript 02b-cellLineLINCS.R /tmp/genes.csv /tmp/drugs.tsv.gz /tmp/cell_line_samples.csv
docker run -v $PWD:/tmp/ cell-line Rscript 02-cellLineDepMap.R /tmp/genes.csv /tmp/cell_line_samples.csv
docker run -v $PWD:/tmp/ cell-line /opt/venv/bin/python 03-curves_and_drug_mapping.py --curSampleFile=/tmp/cell_line_samples.csv
docker run -v $PWD:/tmp cell-line /opt/venv/bin/python 01b-cellLineDrugs_LINCS.py --drugFile drugs.tsv.gz
```

These commands sequentially build the cell line data, samples files, and gene data. 
