# Compass
In-Silico Modeling of Metabolic Heterogeneity using Single-Cell Transcriptomes


## Install

### Installing CPLEX (required for Compass)

Compass depends on IBM CPLEX to run optimization procedures.  The full edition is available for free for academic use.

Install CPLEX [here](https://www.ibm.com/software/commerce/optimization/cplex-optimizer/)

Then follow the insructions [here](https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.0/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html) to install the cplex Python module.


### Installing Compass

The best way is to install directly from this Github repository using the following command:

```bash
pip install git+https://github.com/yoseflab/Compass.git
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
