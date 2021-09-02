Compass Postprocessing
======================

.. toctree::
   :maxdepth: 1
   :caption: Contents:
   :titlesonly:

   notebooks/Demo
   notebooks/Demo-micropools

To demonstrate downstream analysis (postprocessing) of Compass results, we provide a python notebook  `here <https://github.com/YosefLab/Compass/blob/analysis/analysis/Demo.ipynb>`__.
The notebook demonstrates a differential analysis pipeline comparing two groups of cells and replicates figures from our  `paper
<https://doi.org/10.1016/j.cell.2021.05.045>`__ paper analyzing Th17 cell metabolism.


The presentation of the algorithm assumes a single-cell data set.
However, you may choose to group cells together (e.g. via
`metacell <https://github.com/tanaylab/metacell>`__ or
`micropooling <https://github.com/YosefLab/Vision>`__) to reduce
computational overhead. You may also apply Compass to bulk transcriptome
data sets (e.g. bulk RNA-seq or microarray data sets) if there are enough observations (samples) to gain statistical power.

Requirements
************
 - Python
 - Jupyter notebook

All of the required python packages can be installed using a cell at the start of the notebook.

Usage
*****

The example in the notebook uses `metadata <https://github.com/YosefLab/Compass/blob/analysis/analysis/extdata/Th17/cell_metadata.csv>`__ that would come from knowledge of what cells/samples were sequenced 
as well as the `output <https://github.com/YosefLab/Compass/blob/analysis/analysis/extdata/Th17/reactions.tsv>`__ of a Compass run. 

The notebook can be run top down and the main section uses Wilcoxon rank sum to test for differential predicted activity of single reactions between pathogenic Th17 (Th17p) and non-pathogenic Th17 (Th17n) cells.
This differs slightly from the paper's version which tests for predicted differential activity meta-reactions instead, and the code for that is at the end of the notebook.

The notebook can be modified to analyze other datasets by providing a different input of Compass results and relevant cell metadata.
But of course, this example is intended to be a starting point and you may conduct other analyses, replace the Wilcoxon rank sum test or the Cohen's D effect size with your favorite statistical test, and/or add other preprocessing steps depending on how your data looks and what are your research questions.