from __future__ import print_function, division, absolute_import

from .MetabolicModel import (MetabolicModel, Reaction, Compartment,
                     Species, Association, Gene)
from six import string_types
from .geneSymbols import resolve_genes
import libsbml
import re


def load(model_name, sbml_document):

    xml_model = sbml_document.model

    modelx = MetabolicModel(model_name)

    # Add reactions
    for rr in xml_model.getListOfReactions():
        reaction = reaction_from_xml(rr)
        modelx.reactions[reaction.id] = reaction

    # Add compartments
    for cp in xml_model.getListOfCompartments():
        compartment = compartment_from_xml(xml_node=cp)
        modelx.compartments[compartment.id] = compartment

    # Add species
    for met in xml_model.getListOfSpecies():
        species = species_from_xml(xml_node=met)
        modelx.species[species.id] = species

    return modelx


GENE_ASSOCIATION_RE = re.compile("GENE_ASSOCIATION:([^<]*)")


def reaction_from_xml(xml_node):
    """
    Build the reaction from the xml node
    """

    reaction = Reaction()

    reaction.id = xml_node.getId()
    reaction.name = xml_node.getName()
    reaction.subsystem = ""
    reaction.reverse_reaction = None
    reaction.meta = {}

    # Lower and upper bounds
    found_ub = False
    found_lb = False
    for param in xml_node.getKineticLaw().getListOfParameters():
        if param.getId() == 'LOWER_BOUND':
            found_lb = True
            lb = param.getValue()
        elif param.getId() == 'UPPER_BOUND':
            found_ub = True
            ub = param.getValue()

    if not found_lb or not found_ub:
        raise Exception(
            "Error, could not find upper and lower bounds for reaction " +
            reaction.id
        )

    reaction.upper_bound = ub
    reaction.lower_bound = lb

    # Reactants and products

    # Reactants
    reaction.reactants = {}
    for sr in xml_node.getListOfReactants():

        metabolite = sr.getSpecies()
        coefficient = sr.getStoichiometry()

        reaction.reactants.update({
            metabolite: coefficient
        })

    # Products
    reaction.products = {}
    for sr in xml_node.getListOfProducts():

        metabolite = sr.getSpecies()
        coefficient = sr.getStoichiometry()

        reaction.products.update({
            metabolite: coefficient
        })

    # Gene Associations - find the string in the notes
    # Assumes each reaction has a GENE_ASSOCIATION field in its html notes
    gene_association_str = GENE_ASSOCIATION_RE.search(
        xml_node.getNotesString()).group(1).strip()

    if len(gene_association_str) == 0:
        reaction.gene_associations = None
    else:
        reaction.gene_associations = association_from_xml(gene_association_str)

    return reaction


def association_from_xml(rule_str):

    association = _eval_rule_str(rule_str)

    return association


def compartment_from_xml(xml_node):

    compartment = Compartment()

    compartment.id = xml_node.getId()
    compartment.name = xml_node.getName()

    return compartment


def species_from_xml(xml_node):

    species = Species()

    species.id = xml_node.getId()
    species.name = xml_node.getName()
    species.compartment = xml_node.getCompartment()
    species.meta = {}

    return species


# Used to evaluate gene association rules
_TOKEN_RE = re.compile('\(|\)|[^\s\(\)]+')

# Breakdown of RE
"""
    \(       |
    \)       |
    [^\s\(\)]+
"""


def _eval_rule_str(rule):

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

    return _eval_node(elem)


def _eval_node(elem):

    # resolve each node
    # tuple -> association
    # and or or -> operator
    # some other string -> association of type gene
    # end state is list of associations and operators
    resolved = []
    for node in elem:
        if isinstance(node, tuple):
            resolved.append(_eval_node(node))
        elif isinstance(node, string_types):

            if node == "and" or node == "or":
                resolved.append(node)  # must be an operator

            else:
                gene = Gene()

                gene.id = node
                gene.name = node.upper()

                assoc = Association()
                assoc.type = 'gene'
                assoc.gene = gene
                resolved.append(assoc)

        elif isinstance(node, Association):
            resolved.append(node)

        else:
            raise Exception("Unknown Node type: " + str(node))

    if len(resolved) == 1:
        return resolved[0]

    # must be alternating Association and operators
    assert len(resolved) % 2 == 1

    # Look for cases of all | or all & to avoid deep nesting
    found_or = 'or' in resolved
    found_and = 'and' in resolved

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
        # partition on 'or' and recurse

        # Find a middle 'or' to keep trees balanced
        or_indices = [i for i, e in enumerate(resolved) if e == "or"]
        mid = int(len(or_indices)/2)
        i = or_indices[mid]

        left = tuple(resolved[0:i])
        right = tuple(resolved[i + 1:])

        left = _eval_node(left)
        right = _eval_node(right)

        assoc = Association()
        assoc.type = 'or'
        assoc.children = [left, right]
        return assoc

    else:
        raise Exception("Bad node: " + str(resolved))
