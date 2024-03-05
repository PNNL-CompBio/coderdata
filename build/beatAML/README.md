## BeatAML Data generation

This directory builds the data for the BeatAML samples

### Build docker image
```
docker build -f dockerfile.beataml  -t beataml . --build-arg HTTPS_PROXY=$HTTPS_PROXY
```


### Sample generation

```
docker run -v $PWD:/tmp beataml /opt/venv/bin/python GetBeatAML.py --token=$SYNAPSE_AUTH_TOKEN --samples= hcmi_samples.csv

```

### Drug generation

```
docker run -v $PWD:/tmp beataml /opt/venv/bin/python GetBeatAML.py  --token=$SYNAPSE_AUTH_TOKEN


```

### Omics and Experiment Data

```
docker run -v $PWD:/tmp beataml /opt/venv/bin/python GetBeatAML.py  --token=$SYNAPSE_AUTH_TOKEN

```
