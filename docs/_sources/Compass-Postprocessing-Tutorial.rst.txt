Compass Postprocessing
======================

Our `compassR package <https://github.com/YosefLab/compassR>`__
appropriately postprocesses the data and provides an easy, expressive
framework for conducting subsequent analyses based on the
single-cell metabolic reaction consistency matrix produced by the
COMPASS algorithm. It also
includes a suite of expressive utility functions for conducting
statistical analyses building thereupon.

The presentation of the algorithm assumes a single-cell data set.
However, you may choose to group cells together (e.g. via
`metacell <https://github.com/tanaylab/metacell>`__ or
`micropooling <https://github.com/YosefLab/Vision>`__) to reduce
computational overhead. You may also apply COMPASS to bulk transcriptome
data sets (e.g. bulk RNA-seq or microarray data sets) of ample size.

Requirements
************
 - R (3.5.3+)
 - R devtools package

Any devtools version with install_github should work.

CompassR Installation
*********************

1. Make sure you have installed `the ``devtools``
   package <https://github.com/r-lib/devtools>`__ from CRAN.
2. To install the CompassR package from github, run ``devtools::install_github("YosefLab/compassR")``.

You can accomplish both of these steps by running the following R code.

.. code:: r

   # Install devtools from CRAN.
   install.packages("devtools")

   # Install compassR from YosefLab.
   devtools::install_github("YosefLab/compassR")

Usage
*****

In the following tutorial, we’ll explore the Th17 cell data set (`Wagner
et al. <https://www.biorxiv.org/content/10.1101/2020.01.23.912717v1>`__;
`Wang et
al. <https://www.biorxiv.org/content/10.1101/2020.01.23.911966v1>`__)
that ships with the package. It will help you get acquainted with the
basics, while skipping over some of the finer details; if you’re an
advanced user looking for the full documentation, please refer to the
`the wiki <https://github.com/YosefLab/compassR/wiki>`__ instead.

Loading your data
~~~~~~~~~~~~~~~~~

Our first step is to import the libaries we need.

.. code:: r

    library(compassR)
    library(tidyverse)

The next step is to specify a few settings via a ``CompassSettings`` object.

.. code:: r

   compass_settings <- CompassSettings$new(
       user_data_directory = system.file("extdata", "Th17", package = "compassR"),
       cell_id_col_name = "cell_id",
       gene_id_col_name = "HGNC.symbol"
   )

There are 3 important parameters.

-  ``user_data_directory`` is the path to the directory that contains
   the data you want to analyze. This directory should include files
   named ``"cell_metadata.csv"``, ``"reactions.tsv"``, and
   ``"linear_gene_expression_matrix.tsv"``.
-  ``cell_id_col_name`` is the column in ``"cell_metadata.csv"`` that
   uniquely identifies the cells in your data set.
-  And finally, ``gene_id_col_name`` is the column in
   ``"gene_metadata.csv"`` that uniquely identifies the genes you’re
   interested in. Note that you do not have to provide this file unless
   you’re an advanced user. By default, analyses will just use the one
   included in the version of RECON2 that comes with the package – in
   which case you can use ``"HGNC.symbol"`` for human genes or
   ``"MGI.symbol"`` for mouse genes.

Now we can load our data by creating a ``CompassData`` object.

.. code:: r

   compass_data <- CompassData$new(compass_settings)

This line may take a minute to run. Under the hood, it’s postprocessing
the results of the COMPASS algorithm and populating a few tables that
we’ll find useful for our analyses later on:

+-------------------+------+-------------------------------------------+
| Table             | Type | Description                               |
+===================+======+===========================================+
| ``reactio         | Data | Each row is a reaction and each column is |
| n_consistencies`` | f    | a cell. ``reaction_consistencies[i, j]``  |
|                   | rame | is the consitency (or “compatibility”)    |
|                   |      | between reaction ``i`` and cell ``j``.    |
+-------------------+------+-------------------------------------------+
| ``metareactio     | Data | Each row is a metareaction and each       |
| n_consistencies`` | f    | column is a cell.                         |
|                   | rame | ``metareaction_consistencies[i, j]`` is   |
|                   |      | the consistency (or “compatibility”)      |
|                   |      | between metareaction ``i`` and cell       |
|                   |      | ``j``.                                    |
+-------------------+------+-------------------------------------------+
| ``                | Ti   | Each row describes a gene in terms of its |
| metabolic_genes`` | bble | ID and whether it’s a metabolic gene.     |
+-------------------+------+-------------------------------------------+
| ``gene_expres     | Ti   | Each row describes a cell in terms of its |
| sion_statistics`` | bble | ID, total expression, metabolic           |
|                   |      | expression, and metabolic activity.       |
+-------------------+------+-------------------------------------------+
| ``cell_metadata`` | Ti   | The cell metadata from                    |
|                   | bble | ``cell_metadata.csv``. In this example    |
|                   |      | it’s the Th17 cell data from the papers   |
|                   |      | linked above.                             |
+-------------------+------+-------------------------------------------+
| ``gene_metadata`` | Ti   | The gene metadata from the metabolic      |
|                   | bble | model (RECON2, by default).               |
+-------------------+------+-------------------------------------------+
| ``meta            | Ti   | The metabolite metadata from the          |
| bolite_metadata`` | bble | metabolic model (RECON2, by default).     |
+-------------------+------+-------------------------------------------+
| ``re              | Ti   | The reaction metadata from the metabolic  |
| action_metadata`` | bble | model (RECON2, by default).               |
+-------------------+------+-------------------------------------------+
| ``reac            | Ti   | Each row describes a reaction in terms of |
| tion_partitions`` | bble | its ID, undirected ID, direction, and     |
|                   |      | which metareaction (i.e. reaction group)  |
|                   |      | it belongs to.                            |
+-------------------+------+-------------------------------------------+

Note that all the metadata tables’ fields are read as characters, and
must manually be coerced into other data types if desired.

Exploring the statistical analysis suite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we’re ready to start our analysis! We begin by making a
``CompassAnalyzer`` object.

.. code:: r

   compass_analyzer <- CompassAnalyzer$new(compass_settings)

With the ``CompassAnalyzer``, it’s easy to conduct statistical analyses.
Let’s do a Wilcoxon rank-sum test for whether each reaction achieves a
higher consistency among Th17p cells or Th17n cells.

.. code:: r

   group_A_cell_ids <-
       compass_data$cell_metadata %>%
       filter(cell_type == "Th17p") %>%
       pull(cell_id)
   group_B_cell_ids <-
       compass_data$cell_metadata %>%
       filter(cell_type == "Th17n") %>%
       pull(cell_id)
   wilcoxon_results <- compass_analyzer$conduct_wilcoxon_test(
       compass_data$reaction_consistencies,
       group_A_cell_ids,
       group_B_cell_ids,
       for_metareactions = FALSE
   )

We can use functions from the tidyverse to combine the results of our
Wilcoxon test with the data we loaded earlier. Then, with just `a little
``ggplot2`` <ex/>`__, we can even reproduce figures 2(c) and 2(e) from
the papers linked above!
