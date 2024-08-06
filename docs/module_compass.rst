Module-Compass
==============

.. contents:: Contents
   :local:

Module-Compass Settings
*************************

**\-\-select-meta-subsystems** [FILE]
    Invokes Module-Compass and computes Compass scores only for the meta-subsystems listed in the given file.


COMPASS computes reaction scores for each reaction in the entire metabolic network, 
but in many cases users may only be interested in a given collection of subsystems (by which we refer to 
in this documentation as **meta-subsystems**) within the network.

Meta-subsystems are defined in the user-provided textual file. An example input file for Human1 is provided below:

.. code:: bash

    TRANS_META_SUBSYSTEM: Transport reactions
    FA_META_SUBSYSTEM: Fatty acid oxidation; Fatty acid activation (cytosolic); Beta oxidation of even-chain fatty acids (mitochondrial); Fatty acid synthesis; Omega-3 fatty acid metabolism; Omega-6 fatty acid metabolism; Beta oxidation of odd-chain fatty acids (mitochondrial); Fatty acid biosynthesis; Fatty acid biosynthesis (unsaturated); Fatty acid biosynthesis (even-chain); Fatty acid biosynthesis (odd-chain); Beta oxidation of poly-unsaturated fatty acids (mitochondrial); Beta oxidation of unsaturated fatty acids (n-9) (mitochondrial); Fatty acid elongation (even-chain); Beta oxidation of di-unsaturated fatty acids (n-6) (mitochondrial); Beta oxidation of unsaturated fatty acids (n-7) (mitochondrial); Fatty acid metabolism; Beta oxidation of branched-chain fatty acids (mitochondrial); Fatty acid elongation (odd-chain); Fatty acid desaturation (even-chain); Fatty acid desaturation (odd-chain)
    CENTRAL_CARBON_META_SUBSYSTEM: Glycolysis / Gluconeogenesis; Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism
    AA_META_SUBSYSTEM: Arginine and proline metabolism; Phenylalanine, tyrosine and tryptophan biosynthesis; Folate metabolism; Glycine, serine and threonine metabolism; Valine, leucine, and isoleucine metabolism; Tyrosine metabolism; Alanine, aspartate and glutamate metabolism; Cysteine and methionine metabolism; Lysine metabolism; Histidine metabolism; Beta-alanine metabolism; Metabolism of other amino acids; Tryptophan metabolism; Phenylalanine metabolism
    LIPID_META_SUBSYSTEM: Glycerophospholipid metabolism; Glycerolipid metabolism; Steroid metabolism; Glycosphingolipid biosynthesis-lacto and neolacto series; Sphingolipid metabolism; Glycosphingolipid biosynthesis-ganglio series; Cholesterol biosynthesis 1 (Bloch pathway); Eicosanoid metabolism; Cholesterol biosynthesis 3 (Kandustch-Russell pathway); Ether lipid metabolism; Cholesterol metabolism; Glycosphingolipid biosynthesis-globo series; Glycosphingolipid metabolism
    SUGAR_META_SUBSYSTEM: Purine metabolism; Nucleotide metabolism; Pyrimidine metabolism; Starch and sucrose metabolism; Amino sugar and nucleotide sugar metabolism; Fructose and mannose metabolism; Pentose phosphate pathway; Pentose and glucuronate interconversions

Each line in the file defines a meta-subsystem which should be formatted as:

.. code:: bash

    [META_SUBSYSTEM_ID]: [SUBSYSTEM_1]; [SUBSYSTEM_2]; ...

Note that the meta-subsystem IDs should not contain spaces. Also note that the subsystems contained in each meta-subsystem 
should be separated by semicolons.

If you are using Human1, a list of subsystems supported by Module-Compass can be found 
`here <https://github.com/YosefLab/Compass/blob/compass_v2/compass/Resources/Metabolic%20Models/Human1/core_reactions_subsystems.txt>`__.
A list of reactions supported by Module-Compass can be found
`here <https://github.com/YosefLab/Compass/blob/compass_v2/compass/Resources/Metabolic%20Models/Human1/core_reactions_md.csv>`__.

If you are using Mouse1, a list of subsystems supported by Module-Compass can be found 
`here <https://github.com/YosefLab/Compass/blob/compass_v2/compass/Resources/Metabolic%20Models/Mouse1/core_reactions_subsystems.txt>`__.
A list of reactions supported by Module-Compass can be found
`here <https://github.com/YosefLab/Compass/blob/compass_v2/compass/Resources/Metabolic%20Models/Mouse1/core_reactions_md.csv>`__.

If you are using RECON2, a list of subsystems supported by Module-Compass can be found 
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
