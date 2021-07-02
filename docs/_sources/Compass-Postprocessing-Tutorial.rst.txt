Compass Postprocessing
======================

To postprocess the results of the Compass algorithm and analyze them, we have a python notebook  `here <https://github.com/YosefLab/Compass/blob/analysis/analysis/Demo.ipynb>`__ 
which demonstrates a differential analysis pipeline comparing two groups of cells and replicates figures from our paper analyzing Th17 cell metabolism `Wagner
et al. <https://www.biorxiv.org/content/10.1101/2020.01.23.912717v1>`__


The presentation of the algorithm assumes a single-cell data set.
However, you may choose to group cells together (e.g. via
`metacell <https://github.com/tanaylab/metacell>`__ or
`micropooling <https://github.com/YosefLab/Vision>`__) to reduce
computational overhead. You may also apply Compass to bulk transcriptome
data sets (e.g. bulk RNA-seq or microarray data sets) of ample size.

Requirements
************
 - Python
 - Jupyter notebook

All of the required python packages can be installed using a cell at the start of the notebook.

Usage
*****

The example in the notebook uses `metadata <https://github.com/YosefLab/Compass/blob/analysis/analysis/extdata/Th17/cell_metadata.csv>`__ that would come from knowledge of what cells/samples were sequenced 
as well as the `output <https://github.com/YosefLab/Compass/blob/analysis/analysis/extdata/Th17/reactions.tsv>`__ of a Compass run. 

The notebook can be run top down and the main section goes through single-reactions using the Wilcoxon rank-sum test to compare pathogenic Th17p cells to non-pathogenic Th17n cells. 
This differs slightly from the paper's version which uses metareactions instead, and the code for that is at the end of the notebook.

To analyze other datasets you would only need input a different set of Compass results and provide the metadata to identify cells, 
though the figure titles about Th17 cells likely won't apply. 
This example is intended to be a starting point and you may also want to do other analyses, replace the Wilcoxon rank sum test with your favorite statistical test, and/or add other preprocessing steps depending on how your data looks.