Compass User Guide
==================

Welcome to the user guide for Compass, if you haven't already installed Compass visit `the install guide </install>`__. 

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   /install
   /tutorial
   /Compass-Settings
   /Compass-Postprocessing-Tutorial
   /AWS-tutorial

Quick Start
***********

Once Compass has been installed you can run it the algorithm on an input gene expression matrix by opening a command line in the directory with the input file and entering:

.. code:: bash

   compass --data expression.tsv --num-processes 10

If your gene expression matrix is in a sparse matrix market format (.mtx file), you'll also need to provide a tab-delimited file with labels for what genes are in each row

.. code:: bash

   compass --data expression.mtx genes.tsv --num-processes 10

Once Compass has finished running, it will output a matrix ``reactions.tsv`` with reaction penalties where higher scores correspond to a reaction being less likely.

 - For a more detailed example of running Compass, visit `the tutorial </tutorial>`__. 
 - If you want to see all the settings you can customize visit `here </Compass-Settings>`__.
 - For tools to postprocess and analyze the results of Compass, visit `here </Compass-Postprocessing-Tutorial>`__
 - Cloud computing like AWS can help Compass scale because each cell can be handled in parallel, so we provide a `guide <\AWS-tutorial>`__ for running Compass on AWS.