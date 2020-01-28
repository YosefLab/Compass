# Compass
### In-Silico Modeling of Metabolic Heterogeneity using Single-Cell Transcriptomes
Cellular metabolism, a key regulator of immune responses, is difficult to study with current technologies in individual cells. We therefore developed Compass, an in silico approach to infer metabolic status of cells based on transcriptome data, which is readily available even in single cell resolutions.

This page provides an installation guide and a workthrough of a standard analysis. For a detailed description of the algorithm, refer to the biorxiv preprint. (Insert Hyperlink)

## Install

### Installing CPLEX (required for Compass)

Compass depends on IBM CPLEX to run optimization procedures.  The full edition is available for free for academic use.

See our instructions [here](https://github.com/YosefLab/Compass/wiki/Installing-CPLEX-Tutorial) for an in depth tutorial on installing CPLEX

<sub>*Note: Make sure to get the full edition (free for academic use) and not the Community Edition*</sub>

### Installing Compass

The best way is to install directly from this Github repository using the following command:

```bash
pip install git+https://github.com/yoseflab/Compass.git
```

Now to test if everything is installed, simply run:

```bash
compass -h
```
You should see the help text print out if installation was succesful :)

## Running Compass
Broadly speaking, Compass takes in a gene expression matrix scaled by transcripts per million, and outputs a penalty reaction matrix, whereby higher scores correspond to a reaction being **less** likely. 




### Running Compass (Simple)
Input file `expression.txt` should be a tab-delimited text file containing gene expression estimates (TPM) with one row per gene, one column per sample.  Must contain both row and column labels, corresponding to genes and sample IDs/names respectively.

```bash
compass --data expression.txt --num-processes 10
```
Below is an example of the formatting for gene expression (We only show a small portion of the matrix):
![](https://i.imgur.com/cZ4TK9V.jpg)

<sub>*Note: For every individual sample, Compass takes roughly 8 to 24 hours to calculate the reaction penalties (varying by machine). This can be expedited by running more than one process at once.  In addition, Compass saves the results of all samples that it has already processed. Therefore, Compass can also be stopped and restarted after it is done processing a subset of samples.*</sub>

### Running Compass (Advanced settings)
(FILL IN)


### Outputs
When Compass has completed, the outputs for all samples are stored in a tab delimited file `reactions.txt` in the specified output directory (`.` directory when running Compass by default). 

<sub>*Note: While compass is running, it will store partial results for each sample in the `_tmp` directory/*</sub>

