from __future__ import print_function, division, absolute_import

from .MetabolicModel import (MetabolicModel, Reaction, Compartment,
                     Species, Association, Gene)
from ..globals import MODEL_DIR

from six import string_types
import libsbml


def load(model_name, sbml_document, metabolic_model_dir=MODEL_DIR):

    xml_model = sbml_document.model

    # Load Model Parameters
    xml_params = {}
    for pp in xml_model.getListOfParameters():
        key = pp.getId()
        val = pp.getValue()
        xml_params[key] = val

    modelx = MetabolicModel(model_name, metabolic_model_dir=metabolic_model_dir)

    # Add reactions
    for rr in xml_model.getListOfReactions():
        reaction = reaction_from_xml(rr, xml_params)
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


def reaction_from_xml(xml_node, xml_params):
    """
    Build the reaction from the xml node

    Needs xml_params for bounds
    """

    reaction = Reaction()

    reaction.id = xml_node.getId()
    reaction.name = xml_node.getName()
    reaction.subsystem = ""
    reaction.reverse_reaction = None
    reaction.meta = {}

    # Lower and upper bounds

    fbcrr = xml_node.getPlugin('fbc')
    ub = fbcrr.getUpperFluxBound()
    if isinstance(ub, string_types):
        ub = xml_params[ub]

    lb = fbcrr.getLowerFluxBound()
    if isinstance(lb, string_types):
        lb = xml_params[lb]

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

    # Gene Associations
    gpa = fbcrr.getGeneProductAssociation()

    if gpa is None:
        reaction.gene_associations = None
    else:
        reaction.gene_associations = association_from_xml(xml_node=gpa.getAssociation())

    return reaction


def association_from_xml(xml_node):

    association = Association()

    if isinstance(xml_node, libsbml.FbcOr):
        association.type = 'or'
        association.gene = None
        association.children = [
            association_from_xml(xml_node=x)
            for x in xml_node.getListOfAssociations()
        ]

    elif isinstance(xml_node, libsbml.FbcAnd):
        association.type = 'and'
        association.gene = None
        association.children = [
            association_from_xml(xml_node=x)
            for x in xml_node.getListOfAssociations()
        ]

    elif isinstance(xml_node, libsbml.GeneProductRef):
        association.type = 'gene'
        association.children = []

        gp_str = xml_node.getGeneProduct()

        fbmodel = xml_node.getSBMLDocument().model.getPlugin('fbc')
        gp = fbmodel.getGeneProduct(gp_str)

        association.gene = gene_from_xml(xml_node=gp)

    else:
        raise ValueError("Unknown input type: " + type(xml_node))

    return association


def gene_from_xml(xml_node):

    gene = Gene()

    gene.id = xml_node.getId()
    gene.name = xml_node.getName().upper()
    if len(gene.name) == 0:
        gene.name = xml_node.getLabel().upper()
    gene.alt_symbols = []

    return gene


def species_from_xml(xml_node):

    species = Species()

    species.id = xml_node.getId()
    species.name = xml_node.getName()
    species.compartment = xml_node.getCompartment()
    species.formula = xml_node.getPlugin('fbc') \
        .getChemicalFormula()
    species.meta = {}

    return species


def compartment_from_xml(xml_node):

    compartment = Compartment()

    compartment.id = xml_node.getId()
    compartment.name = xml_node.getName()

    return compartment
