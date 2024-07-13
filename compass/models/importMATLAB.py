from __future__ import print_function, division, absolute_import
import os
import json
import pandas as pd
import numpy as np
import re
from ..globals import MODEL_DIR
from .MetabolicModel import Gene, Species, Reaction, Association, MetabolicModel
from six import string_types


def load(model_name, species, metabolic_model_dir=MODEL_DIR):
    """
    model_name: str
        Name of the folder containing the model
    species: str
        Species name.  either 'homo_sapiens' or 'mus_musculus'
    """

    # metabolic_model_dir is 'Resources/Metabolic Models' by default
    # if custom meta-subsystem models are generated, metabolic_model_dir is args['output_dir']/meta_subsystem_models

    # First load Genes
    top_dir = os.path.join(metabolic_model_dir, model_name)
    model_dir = os.path.join(top_dir, 'model')

    with open(os.path.join(model_dir, 'model.genes.json')) as fin:
        gene_ids = json.load(fin)

    with open(os.path.join(top_dir, 'non2uniqueEntrez.json')) as fin:
        gtx = json.load(fin)

    if species == 'homo_sapiens':

        with open(os.path.join(top_dir, 'uniqueHumanGeneSymbol.json')) as fin:
            gene_symbols = json.load(fin)

        alt_symbols = [[] for x in gene_symbols]

    elif species == 'mus_musculus':

        with open(os.path.join(top_dir, 'uniqueMouseGeneSymbol.json')) as fin:
            gene_symbols = json.load(fin)

        with open(os.path.join(top_dir, 'uniqueMouseGeneSymbol_all.json')) as fin:
            alt_symbols = json.load(fin)

    else:
        raise Exception(
            'Unsupported species.  Supported: `homo_sapines`, `mus_musculus`')

    genes = []
    for i, gid in enumerate(gene_ids):

        gene = Gene()

        gene.id = gid
        non_i = gtx[i]-1
        gene.name = gene_symbols[non_i].upper()
        gene.alt_symbols = [x.upper() for x in alt_symbols[non_i]]
        genes.append(gene)

    # Then reactions (evaluate gene rules)
    with open(os.path.join(model_dir, 'model.rxns.json')) as fin:
        rxns = json.load(fin)
    with open(os.path.join(model_dir, 'model.rxnNames.json')) as fin:
        rxnNames = json.load(fin)
    with open(os.path.join(model_dir, 'model.lb.json')) as fin:
        lbs = json.load(fin)
    with open(os.path.join(model_dir, 'model.ub.json')) as fin:
        ubs = json.load(fin)
    with open(os.path.join(model_dir, 'model.subSystems.json')) as fin:
        subSystems = json.load(fin)
    with open(os.path.join(model_dir, 'model.rules.json')) as fin:
        rules = json.load(fin)

    groups = zip(rxns, rxnNames, lbs, ubs, subSystems, rules)
    reactions = []
    for rxn, name, lb, ub, subsystem, rule in groups:

        reaction = Reaction()
        reaction.id = rxn
        reaction.name = name
        reaction.upper_bound = ub
        reaction.lower_bound = lb
        reaction.subsystem = subsystem

        # Eval the rule
        reaction.gene_associations = _eval_rule_str(rule, genes)
        reactions.append(reaction)

    # Other optional reaction files

    # Meta-data
    fname = os.path.join(model_dir, 'rxnMeta.txt')
    if os.path.exists(fname):
        # quoting=3 setting ignores " in file
        rxnMeta = pd.read_csv(fname, sep='\t', index_col=0, quoting=3)

        for i, reaction in enumerate(reactions):
            reaction.meta = _fix_dtypes(
                rxnMeta.loc[reaction.id].to_dict()
            )

    # Then metabolites
    with open(os.path.join(model_dir, 'model.mets.json')) as fin:
        met_ids = json.load(fin)
    with open(os.path.join(model_dir, 'model.metNames.json')) as fin:
        metNames = json.load(fin)
    with open(os.path.join(model_dir, 'model.metFormulas.json')) as fin:
        metFormulas = json.load(fin)

    species = []
    compartment_re = re.compile("\[(.+)\]")
    for met_id, name, formula in zip(met_ids, metNames, metFormulas):

        met = Species()
        met.id = met_id
        met.name = name
        met.formula = formula
        met.compartment = compartment_re.search(met.id).group(1)

        species.append(met)

    # Other optional metabolite files

    # Meta-data
    fname = os.path.join(model_dir, 'metMeta.txt')
    if os.path.exists(fname):
        # quoting=3 setting ignores " in file
        metMeta = pd.read_csv(fname, sep='\t', index_col=0, quoting=3)

        for i, met in enumerate(species):
            met.meta = _fix_dtypes(
                metMeta.loc[met.id].to_dict()
            )

    # Then Smat
    with open(os.path.join(model_dir, 'model.S.json')) as fin:
        Smat = json.load(fin)
        for entry in Smat:
            i = entry[0]
            j = entry[1]
            coef = entry[2]

            species_id = species[i-1].id
            reaction = reactions[j-1]

            if coef > 0:
                reaction.products[species_id] = coef
            else:
                reaction.reactants[species_id] = abs(coef)

    reactions = {r.id: r for r in reactions}
    species = {s.id: s for s in species}
    name = model_name

    model = MetabolicModel(name, metabolic_model_dir=metabolic_model_dir)
    model.reactions = reactions
    model.species = species

    return model

# ----------------------------------------
# Utility functions - used to parse gene association string
# ----------------------------------------

_TOKEN_RE = re.compile('x\(\d+\)|\(|\)|\||&')

# Breakdown of RE
""" x\(\d+\) |
    \(       |
    \)       |
    \|       |
    &
"""


def _eval_rule_str(rule, genes):

    # Return None if there is no gene-product rule
    if len(rule) == 0:
        return None

    # State machine
    token_elem = _TOKEN_RE.findall(rule)

    # Remove parenthesis, replace with tuples
    group_stack = []
    current_group = []
    for term in token_elem:
        if term == '(':
            # Start a new group
            group_stack.append(current_group)
            current_group = []
        elif term == ')':
            # End the group
            prev_group = group_stack.pop()
            prev_group.append(tuple(current_group))
            current_group = prev_group
        else:
            current_group.append(term)

    elem = tuple(current_group)

    return _eval_node(elem, genes)


_ELEM_RE = re.compile('x\((\d+)\)$')


def _eval_node(elem, genes):

    # resolve each node
    # x(2343) -> association of type gene
    # tuple -> association
    # operator -> operator
    # end state is list of associations and operators
    resolved = []
    for node in elem:
        if isinstance(node, tuple):
            resolved.append(_eval_node(node, genes))
        elif isinstance(node, string_types):
            elem_match = _ELEM_RE.match(node)

            if elem_match:
                index = int(elem_match.group(1)) - 1
                gene = genes[index]

                assoc = Association()
                assoc.type = 'gene'
                assoc.gene = gene
                resolved.append(assoc)

            else:
                assert (node == '|' or node == '&')
                resolved.append(node)  # must be an operator
        elif isinstance(node, Association):
            resolved.append(node)

        else:
            raise Exception("Unknown Node type: " + str(node))

    if len(resolved) == 1:
        return resolved[0]

    # must be alternating Association and operators
    assert len(resolved) % 2 == 1

    # Look for cases of all | or all & to avoid deep nesting
    found_or = '|' in resolved
    found_and = '&' in resolved

    if found_or and not found_and:
        nodes = [x for i, x in enumerate(resolved) if i % 2 == 0]

        assoc = Association()
        assoc.type = 'or'
        assoc.children = nodes
        return assoc

    elif not found_or and found_and:
        nodes = [x for i, x in enumerate(resolved) if i % 2 == 0]

        assoc = Association()
        assoc.type = 'and'
        assoc.children = nodes
        return assoc

    elif found_or and found_and:
        # partition on | and recurse

        # Find a middle | to keep trees balanced
        or_indices = [i for i,e in enumerate(resolved) if e == "|"]
        mid = int(len(or_indices)/2)
        i = or_indices[mid]

        left = tuple(resolved[0:i])
        right = tuple(resolved[i + 1:])

        left = _eval_node(left, genes)
        right = _eval_node(right, genes)

        assoc = Association()
        assoc.type = 'or'
        assoc.children = [left, right]
        return assoc

    else:
        raise Exception("Bad node: " + str(resolved))


def _fix_dtypes(input_dict):
    """
    Used to change numpy.int64 and numpy.float64 to python types
    for all values in a dictionary
    """

    out_dict = {}
    for k, v in input_dict.items():
        if isinstance(v, np.int64):
            v = int(v)
        if isinstance(v, np.float64):
            v = float(v)

        out_dict[k] = v

    return out_dict
