# edna-reference-databases
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
