# MFlux
Metabolic Flux-Balance Analysis in Python


## Install

The best way is to install directly from this Github repository using the following command:

```bash
pip install git+https://github.com/yoseflab/mflux.git
```

Note: if you are using `conda` it's best that you first update numpy and pandas through conda by running:

```bash
conda install numpy pandas
```

## Running COMPASS

COMPASS can either be run on a single computer or submitted to a Torque cluster scheduling system.

**Running on a Single Machine**

This will run compass on the current machine using 10 processes.
Input file `expression.txt` should be a tab-delimited text file containing gene expression estimates (TPM) with one row per gene, one column per sample.  Must contain row and column labels.

```bash
compass --data expression.txt --model RECON2_mat --media media1 --num-processes 10
```

**Running through Torque**

Here `queueName` should be the name of the queue to submit the job to.

```bash
compass --data expression.txt --model RECON2_mat --media media1 --torque-queue queueName
```
