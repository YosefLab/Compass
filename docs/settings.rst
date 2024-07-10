Compass Settings
================

.. contents:: Contents
   :local:

Compass allows users to customize various features:

.. code:: bash

   usage: Compass [-h] [--data FILE] [--data-mtx FILE [FILE ...]] [--model MODEL] [--species SPECIES] [--media MEDIA] 
               [--output-dir DIR] [--temp-dir DIR] [--num-processes N] [--torque-queue QUEUE] [--lambda F] [--num-threads N]
               [--and-function FXN] [--select-reactions FILE] [--select-subsystems FILE] [--num-neighbors N] [--symmetric-kernel] 
               [--input-weights FILE] [--penalty-diffusion MODE] [--no-reactions] [--calc-metabolites] [--precache]
               [--input-knn FILE] [--input-knn-distances FILE] [--output-knn FILE] [--latent-space FILE] [--only-penalties]
               [--example-inputs] [--microcluster-size C] [--list-genes FILE] [--list-reactions FILE] [--turbo MIN_SR2]
               [--turbo-increments INC] [--turbo-min-pct MIN_PCT] [--turbo-max-iters MAX_ITERS]


Input Settings
***************

The input gene expression matrix is specified in one of two ways:

**\-\-data** [FILE]
   File with input gene expression data with rows as genes and columns as samples. 
   The input should be a single tab-delimited file with row and column labels:

   .. code:: bash

      --data expression.tsv

**\-\-data-mtx** [--data-mtx FILE [FILE ...]]
   File with input gene expression data with rows as genes and columns as samples in market matrix format (mtx).
   The input must be followed by a tab separated file with rownames corresponding to genes. Optionally that can be followed by column names corresponding to samples.

   .. code:: bash

      --data expression.mtx genes.tsv sample_names.tsv

   If the column names file is omitted the samples will be labelled by index.


To view example inputs, use:

**\-\-example-inputs**
   Flag for Compass to list the directory where example inputs can be found.


Output Settings
****************
   
**\-\-output-dir** [DIR]
   Final directory for final output files (e.g., reactions.tsv). Defaults to ./ (the same directory the command was run from).

**\-\-temp-dir** [DIR]
   Directory to store partial results for completed
   samples in a dataset (used to resume interrupted runs).
   Defaults to ./_tmp.

**\-\-list-genes** [FILE]
   File to output a list of metabolic genes needed for selected metabolic model.
   This is useful if you'd like to subset the input matrix to include only the metabolic genes used by the algorithm
   (gene not included in the list are ignored). This list depends on the ``--species`` argument.
   
**\-\-list-reactions** [FILE]
   File to output a list of reaction id's and their associated subsystem. This is useful if you'd like to compute Compass scores
   for only a subset of the reactions in order to cut in computation times (see below, ``--select-reactions`` and ``--select-subsystems``).

**\-\-select-reactions** [FILE]
   Compute Compass scores only for the reactions listed in the given file. 
   FILE is expected to be textual, with one line per reaction 
   (undirected, namely adding the suffix \"_pos\" or \"_neg\" to a line will create a valid directed reaction id). 
   Unrecognized reactions in FILE are ignored.

**\-\-select-subsystems** [FILE]
   Compute Compass scores only for the subsystems listed in the given file. 
   FILE is expected to be textual, with one line per subsystem.
   Unrecognized subsystems in FILE are ignored.


Metabolic Model Settings
*************************

**\-\-species** [SPECIES]
   Species to use to match gene names to model. Required parameter. Options:

   - homo_sapiens
   - mus_musculus

**\-\-model** [MODEL]
   Metabolic model to use. Options:

   - RECON1_mat 
   - RECON2_mat (default)
   - RECON2.2

**\-\-media** [MEDIA]
   The media to simulate the model with. This is a placeholder for future algorithmic extensions.

**\-\-and-function** [FXN]
   A numeric function which substitutes AND relationships in translation of the GSMM's gene-protein
   associations into reaction penalties Options: 
   
   - min 
   - median
   - mean (default)

**\-\-calc-metabolites**
   Flag to enable calculation and output of
   uptake/secretion scores in addition to reaction scores.

**\-\-no-reactions**
   Flag to disable calculation and output of reaction
   scores and compute only uptake/secretion scores.

Turbo-Compass Settings
***********************

Turbo-Compass is an implementation of Compass that allows for faster runtime at the expense of accuracy. 
If you would like to use Turbo-Compass, please refer to `this section <https://compass-sc.readthedocs.io/en/latest/turbo_compass.html>`__ 
of the documentation.

Penalty Settings
****************

**\-\-penalty-diffusion** [MODE]
   Mode to use in information sharing of reaction penalty values
   between single cells. Options:

   - gaussian 
   - knn (default)

**\-\-lambda** [F]
   Smoothing factor for information sharing between single cells (Default is 0, no information sharing). 
   Lambda should be set between 0 and 1. In the manuscript, where information sharing was appropriate, we used 0.25.
   
   Note there are two common scenarios where there is no need for information sharing and lambda should be set to 0:
   # Running Compass on bulk (i.e., not single cell) RNA
   # Using a cell pooling procedure (`micropools <https://yoseflab.github.io/Compass/micropooling.html>`_, or `metacells <https://tanaylab.github.io/metacell/>`_) and running Compass on the pooled results.
   
.. note::

    If lambda is 0, then the cells are processed independently of each other so you can divide up samples to run them separately and get the same results.

**\-\-num-neighbors** [K]
   Either effective number of neighbors for gaussian
   penalty diffusion or exact number of neighbors for KNN penalty
   diffusion. Default is 30

**\-\-input-weights** [FILE]
   File to input custom weights for averaging of single-cell data.
   The column and row labels should be the same as the names of samples in expression data.

**\-\-symmetric-kernel**
   Flag to enable symmetrizing the TSNE kernel which takes longer


**\-\-input-knn** [FILE]
   File to input a precomputed kNN graph for the samples. 
   File can eiter be a tsv with one row per sample and (k+1) columns. 
   The first column should be sample names, and the next k columns should be indices of the k nearest neighbors (by their order in column 1).

   You can also input the numpy array of values without a column of labels in npy format, but be careful that the order of samples is the same as input data.

**\-\-input-knn-distances** [FILE]
   File to input a precomputed kNN graph for the samples. 
   File can eiter be a tsv with one row per sample and (k+1) columns. 
   The first column should be sample names, and the next k columns should be distances to the k nearest neighbors of that sample.

   You can also input the numpy array of values without a column of labels in npy format, but be careful that the order of samples is the same as input data.
   
**\-\-output-knn** [FILE]
   File to save kNN graph of the samples to.
   File will be a tsv with one row per sample and (k+1) columns. 
   The first column will be sample names, and the next k columns will be indices of the k nearest neighbors (by their order in column 1).

.. note::

   These knn formats are the results from scikit-learn's nearest neighbors algorithm which are then wrapped in a Pandas dataframe.

**\-\-latent-space** [FILE]
   File with latent space representation of the samples on which to do the kNN clustering for information sharing and/or micropooling.
   Should be a tsv with one row per sample and one column per dimension of the latent space.

**\-\-only-penalties**
   Flag for Compass to only compute the reaction penalties for the dataset. This is useful for load splitting when information sharing between cells is needed; only the penalty computation needs to be centrally run, and the subsequent score computations can be split across machines.

Computing Settings
******************

**\-\-num-processes** [N]
   Number of processes for Compass to use, each of which handles a single sample. Must be a positive integer and defaults to the number of processors on machine (using Python's :code:`multiprocessing.cpu_count()`). Ignored
   when submitting job onto a queue

**\-\-num-threads** [N]
   Number of threads to use per sample for solving the flux balance optimization problems. Default is 1. 

.. note::

   It is generally better to increase the number of processes than the number of threads for better performance, unless the number of processes is greater than the number of samples. 
   This is because it is generally better to have multiple optimization problems being solved at once rather than solving a single optimization problem with multiple threads.

**\-\-torque-queue** [QUEUE]
   Name of the torque queue to submit to

**\-\-precache**
   A flag to force Compass to build up the cache for the input selected model and media. This will rebuild the cache even if one already exists.

Microcluster Settings
**********************

.. warning::

   Microclustering is currently deprecated for this version of Compass. We still provide the relevant arguments, but please 
   note that the code is no longer maintained and we do not provide any guarantee on the correctness or validity of the results.
   
   To reduce runtime, we recommend that the user perform pseudobulking on the data, i.e. aggregation of the expression values
   from a group of cells with shared characteristics, such as cells from the same patient, replicate, cell type, etc. 
   As this process is highly dependent on the experiments performed to generate the dataset, 
   we leave this to the discretion of the user.

**\-\-microcluster-size** [C]
   A target number of cells per `microcluster <https://yoseflab.github.io/Compass/micropooling.html>`__. Compass will aggregate similar cells into clusters and compute reaction penalties for the clusters (using the mean of the cluster).

**\-\-microcluster-file** [FILE]
   File where a tsv of microclusters will be output. There will be one column where each entry has the label for what micropool/microcluster the sample is in. Defaults to micropools.tsv in the output directory.

**\-\-microcluster-data-file** [FILE]
   File where a tsv of average gene expression per
   microcluster will be output. Defaults to
   micropooled_data.tsv in the output directory.

.. note::

    When using microclusters, information sharing with lambda > 0 is generally unnecessary because the microclusters already serve the same purpose. If both are enabled, then information will be shared between microclusters as well.

Testing and Hidden Settings
***************************

There are several Compass arguments which are not listed by the parser because they are primarily for testing or for batch jobs.

**\-\-glucose** [G]
   Experimental feature to tweak glucose concentrations in media, default is 1. Higher levels increase glucose availability.

**\-\-test-mode**
   Flag which limits computing scores to the first 100 reactions and first 50 metabolites

**\-\-detailed-perf**
   Flag which enables more performance data collection such as runtimes per reaction per sample.

**\-\-collect** 
   Flag to have Compass collect results. Used for batch jobs

**\-\-config-file** [FILE]
   Setting used for batch jobs

**\-\-penalties-file** [FILE]
   File which allows for specifying a penalties file other than the default one (which is _tmp/penalties.txt.gz)

**\-\-lpmethod** [N]
   Argument to choose the algorithm that Gurobi uses. 
   See the `Gurobi documentation <https://www.gurobi.com/documentation/current/refman/method.html>`__ for more details.
   The default method is automatic selection (-1), although the barrier algorithm (2) and dual simplex (1) also perform well.

**\-\-save-argmaxes**
   Flag to enable saving the argmaxes for computing Compass scores of each reaction. Fun fact: solving the TSP greedily on the argmaxes graph to make full use of the advanced basis setting with the simplex algorithm did not outperform the barrier algorithm in practice.
