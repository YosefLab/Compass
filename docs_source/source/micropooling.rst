Micropooling
============

.. contents:: Contents
   :local:

To help the Compass algorithm scale to larger datasets we use micropooling/microclustering to aggregate samples into clusters which are analyzed with the Compass algorithm which can dramatically lower computational burdens. 
The micropooling algorithm is a reimplementation of micropooling using `VISION <https://www.nature.com/articles/s41467-019-12235-0>`__. 

In short, the algorithm works by first coarsely clustering samples using the Louvain/Leiden algorithm on the K-nearest neighbors (kNN) graph of the gene expression data and then iteratively applying k-means clustering until the clusters are the correct average size.
Each cluster is then treated as a sample where the gene expression is the average of the gene expression for samples in the cluster and the Compass algorithm can be applied to a smaller number of samples. 
If you specify a latent space, the K-nearest neighbor graph will instead be computed using that latent space. You may also specify the kNN graph directly (see in command line arguments).

.. @Brandon: making sure my addition is ok and when you specify kNN directly it is used for micropooling?

Using Micropooling
******************
To enable micropooling using Compass, set the parameter \-\-microcluster-size and the input data will first be clustered into pools before applying the Compass algorithm. 
The micropool assignments will be written to the \-\-microcluster-file in the output directory (by default micropools.tsv) as a tab-delimited table. 
The reactions.tsv file that is output will now be a table of reaction penalties per micropool, analogous to the usual Compass output file but applied to the micropooled data.

To analyze the micropooled results we would instead characterize the micropools. 
* Numeric cell metadata (e.g., library size) can be averaged across cells in the micropool.
* Discrete cell metadata (e.g., belonging to KO or WT) can be decided by majority decision, or if at least 90% (or some other threshold) of the cells in the micropool belong to one of the discrete phenotypes.

 Depending on the research question, it could make sense to micropool discrete phenotypes separately. This will result in micrpools made of only WT or KO cells, for example, but may conceal some of the overlapping cellular programs between the two.


Example of Micropooling
***********************

As in the tutorial, to use microclustering first open a command line in a directory with the gene expression data of interest, say "expression.tsv". You can then run compass on the data with the following command using micropooling:

.. code:: bash

   compass --data expression.tsv --num-processes 10 --microcluster-size 10 --species homo_sapiens

Once the command has finished running, there should be 3 output files in the same directory the command was run: reactions.tsv, micropools.tsv, and micropooled_data.tsv. 
Reactions.tsv contains the penalties per reaction for each micropool, micropools.tsv contains the micropool assignments for each of the original samples, and micropooled_data.tsv contains the average gene expression for each micropool.