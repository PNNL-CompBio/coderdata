## Build test script

To enable future data to be added to CodeData we have built a small
script/files to enable each subsequent docker build to be tested. If
all the commands work then your data is guaranteed to work with the
coderdata framework

### Usage

1. Build your data build package as a docker image, with the
   appropriate scripts (`build_samples.sh`, `build_omics.sh`,
   `build_exp.sh`, `build_drugs.sh`, etc.).
2. From this directory call `python test_docker.py --docker [your
   docker name]` with the following flags:
   - `--samples`
   - `--omics`
   - `--drugs`
   - `--exp`
