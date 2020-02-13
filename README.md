# Compass
### In-Silico Modeling of Metabolic Heterogeneity using Single-Cell Transcriptomes
Cellular metabolism, a key regulator of immune responses, is difficult to study with current technologies in individual cells. We therefore developed Compass, an in silico approach to infer metabolic status of cells based on transcriptome data, which is readily available even in single cell resolutions.

This page provides an installation guide and a workthrough of a standard analysis. For a detailed description of the algorithm, refer to the [biorxiv preprint](https://www.biorxiv.org/content/10.1101/2020.01.23.912717v1?rss=1). 

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
<img src="https://i.imgur.com/VUzHWMr.png" width="646" height="162"></img>


<sub>*Note: For every individual sample, Compass takes roughly 8 to 24 hours to calculate the reaction penalties (varying by machine). This can be expedited by running more than one process at once.  In addition, Compass saves the results of all samples that it has already processed. Therefore, Compass can also be stopped and restarted after it is done processing a subset of samples.*</sub>

### Running Compass (Advanced settings)
Compass also allows users to customize a variaty of settings seen below:
```bash
usage: compass [-h] --data FILE [--model MODEL] [--species SPECIES]
               [--media MEDIA] [--output-dir DIR] [--temp-dir DIR]
               [--torque-queue QUEUE] [--num-processes N] [--lambda F]
               [--num-threads N] [--and-function FXN] [--num-neighbors N]
               [--symmetric-kernel] [--input-weights FILE]
               [--penalty-diffusion MODE] [--no-reactions] [--no-metabolites]
```

See our instructions [here](https://github.com/YosefLab/Compass/wiki/Compass-Advanced-Usage-Tutorial) for an in depth tutorial on using Compass's advanced settings

### Postprocessing

Once Compass has finished running, it is important to apply postprocessing to the data in order to convert reaction penalties(where high values correspond to low likelihood reactions) to reaction scores (Wwhere high values correspond to likely reactions). 

See our instructions [here](https://github.com/YosefLab/Compass/wiki/Compass-Postprocessing-Tutorial) for an in depth tutorial on our postprocessing pipeline.

### Outputs
When Compass has completed, the outputs for all samples are stored in a tab delimited file `reactions.txt` in the specified output directory (`.` directory when running Compass by default). 

Below is an example of the output matrix:
<img src="https://i.imgur.com/pwZbqLT.jpg" width="646" height="143"></img>

<sub>*Note: While compass is running, it will store partial results for each sample in the `_tmp` directory/*</sub>

