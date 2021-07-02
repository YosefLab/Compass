Micropooling
============

.. contents:: Contents
   :local:

To help the Compass algorithm scale to larger datasets we use micropooling/microclustering to aggregate samples into clusters which are analyzed with the Compass algorithm which can dramatically lower computational burdens. 
The micropooling algorithm is a reimplementation of micropooling using `VISION <https://www.nature.com/articles/s41467-019-12235-0>`__. 

In short, the algorithm works by first coarsely clustering samples using the Louvain/Leiden algorithm on the K-nearest neighbors graph of the gene expression data and then iteratively applying k-means clustering until the clusters are the correct average size.
Each cluster is then treated as a sample where the gene expression is the average of the gene expression for samples in the cluster and the Compass algorithm can be applied to a smaller number of samples. 
If you specify a latent space, the K-nearest neighor graph will instead be computed using that latent space.

Using Micropooling
******************
To enable micropooling using Compass, set the parameter \-\-microcluster-size and the input data will first be clustered into pools before applying the Compass algorithm. 
The micropool assignments will be written to the \-\-microcluster-file in the output directory (by default micropools.tsv) as a tab-delimited table. 
The reactions.tsv file that is output will now be a table of reaction penalties per micropool, analagous to the usual Compass output file but applied to the micropooled data.

To analyze the micropooled results we would instead characterize the micropools. For instance when comparing a set of pathogenic cells to non-pathogenic, each micropool would be assessed by what type of cell primarily makes up the micropool.


Example of Micropooling
***********************

As in the tutorial, to use microclustering first open a command line in a directory with the gene expression data of interest, say "expression.tsv". You can then run compass on the data with the following command using micropooling:

.. code:: bash

   compass --data expression.tsv --num-processes 10 --microcluster-size 10 --species homo_sapiens

Once the command has finished running, there should be 3 output files in the same directory the command was run: reactions.tsv, micropools.tsv, and micropooled_data.tsv. 
Reactions.tsv contains the penalties per reaction for each micropool, micropools.tsv contains the micropool assignments for each of the original samples, and micropooled_data.tsv contains the average gene expression for each micropool.