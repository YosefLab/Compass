Tutorial
========

.. contents:: Contents
   :local:

Broadly speaking, Compass takes in a gene expression matrix scaled for library depth 
(e.g., CPM), and outputs a penalty reaction matrix, whereby
higher scores correspond to a reaction being **less** likely.

Running Compass (Simple)
------------------------

The input gene expression matrix can be either a tab-delimited text file (tsv) or a matrix market format (mtx)
containing gene expression estimates (CPM, TPM, or similar scaled units) with one row per gene, one column per sample. 

Tab-delimited files need row and column labels corresponding to genes and sample names. Market matrix formats need a separate tab delimited file of of gene names and optionally a tab delimited file of cell names.

Example input
^^^^^^^^^^^^^

You can find example inputs in tab-delimited format (tsv) and market matrix format (mtx) on this github repo under `compass/Resources/Test-Data <https://github.com/YosefLab/Compass/tree/master/compass/Resources/Test-Data>`__. 

These files will exist locally as well under the Compass install directory which can be found by running:

.. code:: bash

   compass --example-inputs

Running Compass
---------------

Then you can run compass on the data with the following command, which will limit the number of processes to 10:

.. code:: bash

   compass --data expression.tsv --num-processes 10

And to run compass on mtx formatted data use the following:

.. code:: bash

   compass --data-mtx expression.mtx genes.tsv sample_names.tsv --num-processes 10

Though the sample names file can be omitted, in which case the samples will be labelled by index.

Below is an example of the formatting for gene expression (We only show a small portion of the matrix):

.. image:: images/input_ex.png

For the first run of compass on a given model and media there will be overhead building up Compass's cache. 
Compass will automatically build up the cache if it is empty, but you can also manually build up the cache before running compass with:

.. code:: bash

   compass --precache

.. note::
   For every individual sample, Compass takes roughly 30 minutes
   to calculate the reaction penalties (varying by machine). This can
   be expedited by running more than one process at once. In addition,
   Compass saves the results of all samples that it has already processed in the _tmp directory.
   Therefore, Compass can also be stopped and restarted after it is done
   processing a subset of samples so long as the _tmp directory is still there. 

Compass Settings
-----------------

Compass also allows users to customize a variety of settings seen below:

.. code:: bash

   usage: Compass [-h] [--data FILE] [--data-mtx FILE [FILE ...]] [--model MODEL]
               [--species SPECIES] [--media MEDIA] [--output-dir DIR]
               [--temp-dir DIR] [--torque-queue QUEUE] [--num-processes N]
               [--lambda F] [--num-threads N] [--and-function FXN]
               [--select-reactions FILE] [--num-neighbors N]
               [--symmetric-kernel] [--input-weights FILE]
               [--penalty-diffusion MODE] [--no-reactions]
               [--calc-metabolites] [--precache] [--input-knn FILE]
               [--output-knn FILE] [--latent-space FILE] [--only-penalties]
               [--list-genes FILE]


See the instructions
:doc:`here </Compass-Settings>`
for an in depth tutorial on using Compass’s settings

Postprocessing
--------------

Once Compass has finished running, we apply several steps of
postprocessing to the data. Mainly, the postprocessing converts reaction penalties (where
high values correspond to low likelihood reactions) to reaction scores
(where high values correspond to likely reactions).

.. Our `compassR package <https://github.com/YosefLab/compassR>`__
   appropriately post-processes the data and provides an easy, expressive
   framework for conducting subsequent analyses. See :doc:`compass postprocessing tutorial<Compass-Postprocessing-Tutorial>` for more on how to use it.

Outputs
-------

When Compass has completed, the outputs for all samples are stored in a
tab delimited file ``reactions.tsv`` in the specified output directory
(``.`` directory when running Compass by default).

Below is an example of the output matrix:

.. image:: images/output_ex.png

To get more context on what the reaction identifiers are, you can visit `virtual metabolic human <https://www.vmh.life/#home>`__ or `the resources directory of compass <https://github.com/YosefLab/Compass/tree/master/compass/Resources/Recon2_export>`__ where there are several csv's which include information on the reactions in Recon2.

\ *Note: While compass is running, it will store partial results for
each sample in the _tmp directory/ (or the directory following \-\-temp\-dir)*\ 
