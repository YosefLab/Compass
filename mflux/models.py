"""
For working with metabolic models
"""
from __future__ import print_function, division
import os
import libsbml
from .globals import RESOURCE_DIR

MODEL_DIR = os.path.join(RESOURCE_DIR, 'Metabolic Models')


def load_metabolic_model(file_name):
    """
    Loads the metabolic model from `file_name`, returning a Model object
    """

    if file_name.endswith('.xml'):
        return load_metabolic_model_xml(file_name)
    else:
        raise NotImplementedError("Can only handle xml currently")


def load_metabolic_model_xml(file_name):

    full_path = os.path.join(MODEL_DIR, file_name)
    sbmlDocument = libsbml.readSBMLFromFile(full_path)
    xml_model = sbmlDocument.model

    # Load Model Parameters
    xml_params = {}
    for pp in xml_model.getListOfParameters():
        key = pp.getId()
        val = pp.getValue()
        xml_params[key] = val

    modelx = MetabolicModel()

    # Add reactions
    for rr in xml_model.getListOfReactions():
        reaction = Reaction(xml_node=rr, xml_params=xml_params)
        modelx.reactions[reaction.id] = reaction

    # Add compartments
    for cp in xml_model.getListOfCompartments():
        compartment = Compartment(xml_node=cp)
        modelx.compartments[compartment.id] = compartment

    # Add species
    for met in xml_model.getListOfSpecies():
        species = Species(xml_node=met, model=modelx)
        modelx.species[species.id] = species

    return modelx


class MetabolicModel(object):

    def __init__(self):
        self.reactions = {}
        self.species = {}
        self.compartments = {}
        self.objectives = {}

    def getReactions(self):
        """
        Returns a list of reaction id's in the MetabolicModel
        """
        return list(self.reactions.keys())

    def getReactionBounds(self):
        """
        Returns two dicts, each mapping reaction id -> bound

        Returns:
            lower_bounds
            upper_bounds
        """

        lower_bounds = {x.id: x.lower_bound for x in self.reactions.values()}
        upper_bounds = {x.id: x.upper_bound for x in self.reactions.values()}

        return lower_bounds, upper_bounds

    def getReactionExpression(self, expression):
        # type: (pandas.Series) -> dict
        """
        Evaluates a score for every reaction, using the expression data.
        This is used for constraint/penalty generation.
        If a score cannot be computed, NaN is used to fill

        Returns a dict:
            key: reaction id (str)
            val: reaction score (float)
        """

        score_dict = {}

        for r_id, reaction in self.reactions.items():
            score = reaction.eval_expression(expression)
            score_dict[r_id] = score

        return score_dict

    def limitUptakeReactions(self, limit):
        """
        Limits the rate of metabolite exchange.

        Applies the limit of `limit` where the limit is non-zero
        Directionality is preserved

        Exchange reactions are defined as any reaction in which
        metabolites are produced from nothing

        Does not return anything - modifies objects in place
        """

        exchanges = {r_id: r for r_id, r in self.reactions.items()
                     if r.is_exchange}

        # Make sure all exchange reactions have reactants but no products
        # Instead of the opposite
        for reaction in exchanges.values():
            if len(reaction.products) > 0 and len(reaction.reactants) == 0:
                reaction.invert()

        # Negative flux now means that metabolites are created
        # Limit this rate
        for reaction in exchanges.values():

            if reaction.lower_bound < -1 * limit:
                reaction.lower_bound = -1 * limit

            if reaction.upper_bound < -1 * limit:
                reaction.upper_bound = -1 * limit

    def getSMAT(self):
        """
        Returns a sparse form of the s-matrix

        result is a dict
            key: metabolite (species) id
            value: list of 2-tuples (reaction_id, coefficient)

        coefficient is positive if metabolite is produced in the reaction,
            negative if consumed
        """
        s_mat = {}
        for reaction_id, rr in self.reactions.items():

            # reactants
            for metabolite, coefficient in rr.reactants.items():

                if metabolite not in s_mat:
                    s_mat[metabolite] = []

                s_mat[metabolite].append((reaction_id, coefficient * -1))

            # products
            for metabolite, coefficient in rr.products.items():

                if metabolite not in s_mat:
                    s_mat[metabolite] = []

                s_mat[metabolite].append((reaction_id, coefficient))

        return s_mat

    def make_unidirectional(self):
        """
        Splits every reaction into a positive and negative counterpart
        """
        uni_reactions = {}
        for reaction in self.reactions.values():

            # Copy the pos_reaction
            pos_reaction = Reaction(from_reaction=reaction)

            # Copy the negative reaction and Invert
            neg_reaction = Reaction(from_reaction=reaction)
            neg_reaction.invert()

            # Adjust bounds for both
            # Examples
            # Positive: Original -> Clipped
            #             0:10   -> 0:10
            #           -10:0    -> 0:0
            #             5:7    -> 5:7
            #            -9:-5   -> 0:0
            #
            # Negative: Original -> Flipped -> Clipped
            #             0:10   -> -10:0   -> 0:0
            #           -10:0    ->   0:10  -> 0:10
            #             5:7    ->  -7:-5  -> 0:0
            #            -9:-5   ->   5:9   -> 5:9

            if pos_reaction.upper_bound < 0:
                pos_reaction.upper_bound = 0

            if pos_reaction.lower_bound < 0:
                pos_reaction.lower_bound = 0

            if neg_reaction.upper_bound < 0:
                neg_reaction.upper_bound = 0

            if neg_reaction.lower_bound < 0:
                neg_reaction.lower_bound = 0

            pos_reaction.id = pos_reaction.id + "_pos"
            neg_reaction.id = neg_reaction.id + "_neg"

            uni_reactions[pos_reaction.id] = pos_reaction
            uni_reactions[neg_reaction.id] = neg_reaction

        self.reactions = uni_reactions

        # Remove reactions that can't carry flux
        self.reactions = {r_id: r for r_id, r in self.reactions.items()
                          if r.upper_bound > 0}


class Reaction(object):
    # Bounds
    # Reactants (and coefficients)
    # Products (and coefficients)
    # Gene Associations

    def __init__(self, xml_node=None, xml_params=None,
                 from_reaction=None):

        if xml_node is not None and xml_params is not None:
            Reaction.__init__xml(self, xml_node, xml_params)

        elif from_reaction is not None:  # Copy constructor

            self.id = from_reaction.id
            self.upper_bound = from_reaction.upper_bound
            self.lower_bound = from_reaction.lower_bound
            self.reactants = from_reaction.reactants.copy()
            self.products = from_reaction.products.copy()
            self.gene_associations = from_reaction.gene_associations

        else:  # Placeholders

            self.id = ""
            self.upper_bound = float('nan')
            self.lower_bound = float('nan')
            self.reactants = {}
            self.products = {}
            self.gene_associations = None

    def __init__xml(self, xml_node, xml_params):
        """
        Build the reaction from the xml node

        Needs xml_params for bounds
        """
        self.id = xml_node.getId()

        # Lower and upper bounds

        fbcrr = xml_node.getPlugin('fbc')
        ub = fbcrr.getUpperFluxBound()
        if isinstance(ub, str):
            ub = xml_params[ub]

        lb = fbcrr.getLowerFluxBound()
        if isinstance(lb, str):
            lb = xml_params[lb]

        self.upper_bound = ub
        self.lower_bound = lb

        # Reactants and products

        # Reactants
        self.reactants = {}
        for sr in xml_node.getListOfReactants():

            metabolite = sr.getSpecies()
            coefficient = sr.getStoichiometry()

            self.reactants.update({
                metabolite: coefficient
            })

        # Products
        self.products = {}
        for sr in xml_node.getListOfProducts():

            metabolite = sr.getSpecies()
            coefficient = sr.getStoichiometry()

            self.products.update({
                metabolite: coefficient
            })

        # Gene Associations
        gpa = fbcrr.getGeneProductAssociation()

        if gpa is None:
            self.gene_associations = None
        else:
            self.gene_associations = Association(xml_node=gpa.getAssociation())

    def eval_expression(self, expression):
        """
        Evalutes a reaction expression.
        Returns 'nan' if there is no gene association
        """

        if self.gene_associations is None:
            return float('nan')
        else:
            return self.gene_associations.eval_expression(expression)

    def invert(self):
        """
        Inverts the directionality of the reaction in-place
        """

        products = self.products
        reactants = self.reactants

        self.reactants = products
        self.products = reactants

        lower = self.lower_bound
        upper = self.upper_bound

        self.upper_bound = lower * -1
        self.lower_bound = upper * -1

    @property
    def is_exchange(self):

        if len(self.products) == 0 and len(self.reactants) > 0:
            return True
        elif len(self.reactants) == 0 and len(self.products) > 0:
            return True
        else:
            return False


class Species(object):

    def __init__(self, model, xml_node=None):
        if xml_node is not None:
            Species.__init__xml(self, model, xml_node)
        else:
            self.id = ""
            self.compartment = ""

    def __init__xml(self, model, xml_node):

        self.id = xml_node.getId()
        self.compartment = model.compartments[
            xml_node.getCompartment()]


class Compartment(object):

    def __init__(self, xml_node=None):
        if xml_node is not None:
            Compartment.__init__xml(self, xml_node)
        else:
            self.id = ""
            self.name = ""

    def __init__xml(self, xml_node):

        self.id = xml_node.getId()
        self.name = xml_node.getName()


class Association(object):

    def __init__(self, xml_node=None):
        if xml_node is not None:
            Association.__init__xml(self, xml_node)
        else:
            self.type = ''
            self.gene = None
            self.children = []

    def __init__xml(self, xml_node):

        if isinstance(xml_node, libsbml.FbcOr):
            self.type = 'or'
            self.gene = None
            self.children = [Association(xml_node=x)
                             for x in xml_node.getListOfAssociations()]

        elif isinstance(xml_node, libsbml.FbcAnd):
            self.type = 'and'
            self.gene = None
            self.children = [Association(xml_node=x)
                             for x in xml_node.getListOfAssociations()]

        elif isinstance(xml_node, libsbml.GeneProductRef):
            self.type = 'gene'
            self.children = []

            gp_str = xml_node.getGeneProduct()

            fbmodel = xml_node.getSBMLDocument().model.getPlugin('fbc')
            gp = fbmodel.getGeneProduct(gp_str)

            self.gene = Gene(xml_node=gp)

        else:
            raise ValueError("Unknown input type: " + type(xml_node))

    def eval_expression(self, expression):
        """
        Matches genes with scores
        Combines 'or' associations with 'sum'
        Combines 'and' associations with 'min'
        """

        if self.type == 'and':
            return min([x.eval_expression(expression) for x in self.children])

        elif self.type == 'or':
            return sum([x.eval_expression(expression) for x in self.children])

        elif self.type == 'gene':

            try:
                return expression[self.gene.name]
            except KeyError:
                return 0.0

        else:
            raise Exception("Unknown Association Type")


class Gene(object):

    def __init__(self, xml_node=None):
        if xml_node is not None:
            Gene.__init__xml(self, xml_node)
        else:
            self.id = ''
            self.name = ''

    def __init__xml(self, xml_node):

        self.id = xml_node.getId()
        self.name = xml_node.getName().upper()
