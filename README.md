# edna-reference-databases
## Releases
### 20231101

:link: [20231101.zip](https://edna-reference-databases.s3.us-east-1.amazonaws.com/20231101.zip)

## Dependencies

```bash
conda create -n crabs
conda activate crabs
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda crabs
conda install -c bioconda rdptools
conda install -c bioconda seqtk
```

## Run

```bash
python -m ednarefdb
```
