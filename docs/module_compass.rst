Module-Compass
==============

.. contents:: Contents
   :local:

Module-Compass Settings
*************************

COMPASS computes reaction scores for each reaction in the entire metabolic network, 
but in many cases users may only be interested in a given collection of subsystems (by which we refer to 
in this documentation as **meta-subsystems**) within the network. To invoke Module-Compass, 
run COMPASS with the flag ``--select-meta-subsystems [PATH_TO_META_SUBSYSTEMS_FILE]``.

Meta-subsystems are defined in the user-provided textual file. An example input file is provided below:

.. code:: bash

    Extracellular transport (TRANS_EXT_META_PTHWY): Transport, extracellular
    Fatty acid metabolism (FA_META_PTHWY): Fatty acid oxidation; Fatty acid synthesis
    Glycolysis (GLY_META_PTHWY): Glycolysis/gluconeogenesis
    Citrid acid cycle (TCA_META_PTHWY): Citric acid cycle

Each line in the file defines a meta-subsystem which should be formatted as:

.. code:: bash

    [META_SUBSYSTEM_NAME] ([META_SUBSYSTEM_ID]): [SUBSYSTEM_1]; [SUBSYSTEM_2]; ...

Note that the meta-subsystem IDs should not contain spaces and should be formatted in parentheses. Also note that 
the subsystems contained in each meta-subsystem should be separated by semicolons. 

A list of subsystems supported by Module-Compass can be found 
`here <https://github.com/YosefLab/Compass/blob/compass_v2/compass/Resources/Metabolic%20Models/RECON2_mat/model/core_reactions_subsystems.txt>`__.

A list of reactions supported by Module-Compass can be found
`here <https://github.com/YosefLab/Compass/blob/compass_v2/compass/Resources/Metabolic%20Models/RECON2_mat/model/core_reactions_md.csv>`__.


Module-Compass Explained
**************************

To support running COMPASS on a subset of the network, we hereby provide Module-Compass, an algorithm
that treats user-defined meta-subsystems as individual networks but still provides the necessary context 
to ensure that this simplification of the network is reasonable.

More specifically, each meta-subsystem is defined as a collection of subsystems that the user chooses. 
We construct the network for each meta-subsystem by first taking all reactions associated with the subsystems within 
the given meta-subsystems, then taking their one-step neighbor reactions, and finally adding exchange reactions for 
all metabolites associated with the collection of reactions we constructed above. The neighboring reactions provides 
necessary context to the meta-subsystem to ensure a reasonable output, while the exchange reactions enables us to 
treat the constructed network as one that is able to interface with the surrounding environment 
and uptake/secrete metabolites.

One thing to note is that the ``--select-reactions`` and ``--select-subsystems`` flags also allow users to specify 
a list of reactions and subsystems to compute COMPASS on. However, Module-Compass differs fundamentally from these 
methods. ``--select-reactions`` and ``--select-subsystems`` both operate on the entire metabolic network, meaning that the 
computation of each reaction score is still constrained by all reactions in the entire metabolic network. This requires 
the linear solver to optimize for all variables (reactions). On the other hand, Module-Compass treats meta-subsystems as 
standalone networks, drastically reducing the number of variables in the linear optimization problem and therefore 
resulting in significant speedups.
