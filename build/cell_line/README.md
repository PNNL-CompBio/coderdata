## Building cell line data
The cell line data is the first to be built, and requires the following commands.

```
docker build -f ../dockerfile.cell-line -t cell-line ../
docker run cell-line Rscript 00-buildGeneFile.R
docker run cell-linne Rscript 01-cellLineSamples.R
docker run cell-line Rscript 02-cellLineDepMap.R genes.csv cell_line_samples.csv
docker run cell-line Rscript 02a-cellLineSanger.R genes.csv cell_line_samples.csv
docker run cell-line python 03-curves_and_drug_mapping.py --curSampleFile=cell_line_samples.csv
```

These commands sequentially build the cell line data, samples files, and gene data. 
