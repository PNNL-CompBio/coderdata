## Building CPTAC data

To build the CPTAC data execute the following commands..

```
docker build -f dockerfile.cptac  -t cptac . --build-arg HTTPS_PROXY=$HTTPS_PROXY

docker run -v $PWD:/tmp cptac --geneFile=genes.csv --prevSampleFile=/tmp/cell_line_samples.csv 
docker run -v $PWD:/tmp cptac --geneFile=genes.csv --curSampleFile=/tmp/cptac_samples.csv 


```
