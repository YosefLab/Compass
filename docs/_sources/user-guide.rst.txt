Compass User Guide
==================

Welcome to the user guide for Compass, if you haven't already installed Compass visit :doc:`the install guide <install>`. 

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

If your gene expression matrix is in a matrix market format (.mtx file) you can use the ``--data-mtx`` option. You'll also need to provide a tab-delimited file with labels for what genes are in each row and optionally add in a sample names file:

.. code:: bash

   compass --data-mtx expression.mtx genes.tsv sample_names.tsv --num-processes 10

Once Compass has finished running, it will output a matrix ``reactions.tsv`` with reaction penalties where higher scores correspond to a reaction being less likely.


 - For a more detailed example of running Compass, visit :doc:`the tutorial <tutorial>`. 
 - If you want to see all the settings you can customize visit :doc:`here <Compass-Settings>`.
 - For tools to postprocess and analyze the results of Compass, visit :doc:`here </Compass-Postprocessing-Tutorial>`
 - Cloud computing like AWS can help Compass scale because each cell can be handled in parallel, so we provide a :doc:`guide <AWS-tutorial>` for running Compass on AWS.