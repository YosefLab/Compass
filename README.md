# Compass
### In-Silico Modeling of Metabolic Heterogeneity using Single-Cell Transcriptomes
Cellular metabolism is a major regulator of immune response, but it is not easy to study the metabolic status of an individual immune cell using current technologies. This motivated us to develop an in silico approach to infer metabolic status of an immune cell by using single-cell transcriptomes. Here, we present COMPASS, an algorithm to characterize the metabolic landscape of single cells based on single-cell RNA-Seq profiles and flux balance analysis.

For specific details regarding our method, please see the following manuscript (Insert Hyperlink)

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

## Running COMPASS
Broadly speaking, Compass takes in a gene expression matrix scaled by transcripts per million, and outputs a penalty reaction matrix, whereby higher scores correspond to a reaction being **less** likely. 




### Running Compass (Simple)
Input file `expression.txt` should be a tab-delimited text file containing gene expression estimates (TPM) with one row per gene, one column per sample.  Must contain both row and column labels, corresponding to genes and sample IDs/names respectively.

```bash
compass --data expression.txt --num-processes 10
```

<sub>*Note: For every individual sample, Compass takes roughly 8 to 24 hours to calculate the reaction penalties (varying by machine). This can be expedited by running more than one process at once.  In addition, Compass saves the results of all samples that it has already processed. Therefore, Compass can also be stopped and restarted after it is done processing a subset of samples.*</sub>

### Running Compass (Advanced settings)
(FILL IN)


### Outputs
When COMPASS has completed, the outputs for all samples are stored in a tab delimited file 'reactions.txt' in the specified output directory ('.' directory when running COMPASS by default). 

<sub>*Note: While compass is running, it will store partial results for each sample in the '_tmp' directory/*</sub>

