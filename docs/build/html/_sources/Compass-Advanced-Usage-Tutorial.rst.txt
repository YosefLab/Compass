Advanced Usage Tutorial
=======================

Compass allows users to customize various features when using compass.

.. code:: bash

   usage: compass [-h] --data FILE [--model MODEL] [--species SPECIES]
                  [--media MEDIA] [--output-dir DIR] [--temp-dir DIR]
                  [--torque-queue QUEUE] [--num-processes N] [--lambda F]
                  [--num-threads N] [--and-function FXN] [--num-neighbors N]
                  [--symmetric-kernel] [--input-weights FILE]
                  [--penalty-diffusion MODE] [--no-reactions] [--no-metabolites]

Below we describe the features in more detail:

-  **and-function**: Which function used to aggregate AND
   associations(min, median, mean). Default min
-  **input-weights**:
-  **lambda**: Smoothing factor for single-cell data. Should be set
   between 0 and 1
-  **media**:
-  **model**: Metabolic model to use (RECON1_mat, RECON2_mat, RECON2.2).
   Default is RECON2_mat
-  **no-metabolites**: Flag to disable calculation and output of
   uptake/secretion scores in addition to reaction scores.
-  **no-reactions**: Flag to disable calculation and output of reaction
   scores in addition to uptake/secretion scores.
-  **num-neighbors**: Either effective number of neighbors for gaussian
   penalty diffusion or exact number of neighbors for KNN penalty
   diffusion. Default is 30
-  **num-processes**: Number of processes for Compass to use. Ignored
   when submitting job onto a queue
-  **num-threads**: Number of threads to use per sample. Default is 1. Generally scaling the number of processes results is better than the number of threads.
-  **output-dir**: Final directory to output concatenated reactions.txt
   file
-  **penalty-diffusion**: Mode to use to share reaction penalty values
   between single cells (gaussian, knn). Default gaussian
-  **species**: Species to use to match genes to model (homo_sapiens,
   mus_musculus). Default is homo_sapiens
-  **symmetric-kernel**:
-  **temp-dir**: Directory to store partial results for completed
   samples in a dataset
-  **torque-queue**: Name of the torque queue to submit to
