Compass Settings
================

.. contents:: Contents
   :local:

Compass allows users to customize various features:

.. code:: bash

   usage: Compass [-h] [--data FILE] [--data-mtx FILE [FILE ...]] [--model MODEL] [--species SPECIES] [--media MEDIA] [--output-dir DIR]
               [--temp-dir DIR] [--torque-queue QUEUE] [--num-processes N] [--lambda F] [--num-threads N] [--and-function FXN]
               [--select-reactions FILE] [--select-subsystems FILE] [--num-neighbors N] [--symmetric-kernel] [--input-weights FILE]
               [--penalty-diffusion MODE] [--no-reactions] [--calc-metabolites] [--precache] [--input-knn FILE] [--output-knn FILE]
               [--latent-space FILE] [--only-penalties] [--example-inputs] [--microcluster-size C] [--list-genes FILE] [--list-reactions FILE]


Below we describe the features in more detail:

Input and Output settings
-------------------------

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

**\-\-output-dir** [DIR]
   Final directory to output concatenated reactions.txt file

**\-\-temp-dir** [DIR]
   Directory to store partial results for completed
   samples in a dataset

**\-\-example-inputs**
   Flag for Compass to list the directory where example inputs can be found.

Metabolic Model Settings
------------------------

**\-\-model** [MODEL]
   Metabolic model to use. Options:

   - RECON1_mat 
   - RECON2_mat (default)
   - RECON2.2

**\-\-media** [MEDIA]
   The media to simulate the model with.

**\-\-species** [SPECIES]
   Species to use to match genes to model. Options:

   - homo_sapiens (default)
   - mus_musculus

**\-\-and-function** [FXN]
   Which function used to aggregate and
   associations. Options: 
   
   - min 
   - median
   - mean (default)

**\-\-calc-metabolites**
   Flag to enable calculation and output of
   uptake/secretion scores in addition to reaction scores.

**\-\-no-reactions**
   Flag to disable calculation and output of reaction
   scores in addition to uptake/secretion scores.

**\-\-list-genes** [FILE]
   File to output a list of metabolic genes needed for selected metabolic model.
   
**\-\-list-reactions** [FILE]
   File to output a list of reaction id's and their associated subsystem for selected metabolic model as a json file.

**\-\-select-reactions** [FILE]
   Compute compass scores only for the reactions listed in the given file. 
   FILE is expected to be textual, with one line per reaction 
   (undirected, namely adding the suffix \"_pos\" or \"_neg\" to a line will create a valid directed reaction id). 
   Unrecognized reactions in FILE are ignored.

**\-\-select-subsystems** [FILE]
   Compute compass scores only for the subsystems listed in the given file. 
   FILE is expected to be textual, with one line per subsystem.
   Unrecognized subsystem in FILE are ignored.


Penalty Settings
----------------

**\-\-penalty-diffusion** [MODE]
   Mode to use to share reaction penalty values
   between single cells. Options:

   - gaussian 
   - knn (default)

**\-\-lambda** [F]
   Smoothing factor for single-cell data. Should be set between 0 and 1

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
   File must be a tsv with one row per sample and (k+1) columns. 
   The first column should be sample names, and the next k columns should be indices of the k nearest neighbors (by their order in column 1).

**\-\-output-knn** [FILE]
   File to save kNN graph of the samples to.
   File will be a tsv with one row per sample and (k+1) columns. 
   The first column will be sample names, and the next k columns will be indices of the k nearest neighbors (by their order in column 1).

**\-\-latent-space** [FILE]
   File with latent space representation of the samples on which to do the kNN clustering.
   Should be a tsv with one row per sample and one column per dimension of the latent space.

**\-\-only-penalties**
   Flag for Compass to only compute the reaction penalties for the dataset.

Computing Settings
------------------

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
   A flag to force compass to build up the cache for the input selected model and media. This will rebuild the cache even if one already exists.

**\-\-microcluster-size** [C]
   A target number of cells per microcluster. Compass will aggregate similar cells into clusters and compute reaction penalties for the clusters (using the mean of the cluster).

.. note::
    When using microclusters, information sharing with lambda > 0 is generally unneccesary because the microclusters already serve the same purpose. If both are enabled, then information will be shared between microclusters as well.

Testing and Hidden Settings
---------------------------
There are several Compass arguments which are not listed by the parser because they are primarily for testing or for batch jobs.

**\-\-glucose** [G]
   Experimental feature to tweak glucose concentrations in media, default is 1. Higher levels increase glucose availability.

**\-\-test-mode**
   Flag which limits computing scores to the first 100 reactions and first 50 metabolites

**\-\-detailed-perf**
   Flag which enables more performance data collection such as runtimes per reaction per sample.

**\-\-collect** 
   Flag to have compass collect results. Used for batch jobs

**\-\-config-file** [FILE]
   Setting used for batch jobs

**\-\-penalties-file** [FILE]
   File which allows for specifying a penalties file other than the default one (which is _tmp/penalties.txt.gz)

**\-\-lpmethod** [N]
   Argument to choose the algorithm CPLEX uses. 
   See `Cplex documentation for more details <https://www.ibm.com/support/knowledgecenter/SSSA5P_20.1.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/LPMETHOD.html>`__. 
   Through testing the barrier algorithm (4) is fastest and therefore default, with automatic selection (0) or dual simplex (2) also performing well.

**\-\-advance** [N]
   Argument to choose the setting for Cplex's advanced basis setting.
   See `Cplex documentaton for more details <https://www.ibm.com/support/knowledgecenter/SSSA5P_20.1.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/AdvInd.html>`__.
   Defaults to 2 as best runtime was found using that for tests.

**\-\-save-argmaxes**
   Flag to enable saving the argmaxes for computing Compass scores of each reaction. Fun fact: solving the TSP greedily on the argmaxes graph to make full use of the advanced basis setting with the simplex algorithm did not outperform the barrier algorithm in practice.




