## HCMI Data

Here we will store the scripts required to process the data from the [Human Cancer Models Initiative](https://ocg.cancer.gov/programs/HCMI)

## Build Docker

```
docker build -f build/docker/Dockerfile.pancpdo -t pancpdo . 
```

## Run build command
```
python build_dataset.py --dataset pancpdo --build
```

