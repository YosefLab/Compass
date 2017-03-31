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

    def getReactionScores(self, expression):
        # type: (pandas.Series) -> dict
        """
        Evaluates a score for every reaction, using the expression data.
        This is used for constraint generation.
        If a score cannot be computed, NaN is used to fill

        Returns a dict:
            key: reaction id (str)
            val: reaction score (float)
        """

        score_dict = {}

        for r_id, reaction in self.reactions.items():
            score = reaction.eval_score(expression)
            score_dict[r_id] = score

        return score_dict

    def limitUptakeReactions(self, limit):
        """
        Limits the rate of metabolite exchange.

        Applies the limit of `limit` where the limit is non-zero
        Directionality is preserved

        Exchange reactions are defined as any reaction in which
        one metabolite is in the extracellular compartment and another is not

        Does not return anything - modifies objects in place
        """

        ex_id = self._getExtracellularCompartmentId()

        for rr in self.reactions.values():

            species = list(rr.reactants.keys())
            species = species + list(rr.products.keys())

            species = [self.species[x] for x in species]

            compartments = [x.compartment for x in species]

            # Count # of ex_id
            n = 0
            for c in compartments:
                if c.id == ex_id:
                    n += 1

            if n > 0 and n < len(compartments):  # This is an exchange reaction

                old_lb = rr.lower_bound
                if old_lb < 0:
                    rr.lower_bound = -1 * limit
                elif old_lb > 0:
                    rr.lower_bound = limit
                # Do not modify if limit was equal to 0

                old_ub = rr.upper_bound
                if old_ub < 0:
                    rr.upper_bound = -1 * limit
                elif old_ub > 0:
                    rr.upper_bound = limit
                # Do not modify if limit was equal to 0

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

    def getExtracellularMetabolites(self):
        """
        Returns a set of the IDs of extracellular metabolites
        """

        ex_id = self._getExtracellularCompartmentId()

        ex_metabs = set()
        for met_id, met in self.species.items():
            if met.compartment.id == ex_id:
                ex_metabs.add(met_id)

        return ex_metabs

    def _getExtracellularCompartmentId(self):
        """
        Identifies the extracellular compartment and returns its ID
        """

        ex_id = []

        for cp in self.compartments.values():
            if 'extracellular' in cp.name:
                ex_id.append(cp.id)

        if len(ex_id) != 1:
            raise Exception("Can't Determine Extracellular "
                            "Compartment in model")

        ex_id = ex_id[0]

        return ex_id


class Reaction(object):
    # Bounds
    # Reactants (and coefficients)
    # Products (and coefficients)
    # Gene Associations

    def __init__(self, xml_node=None, xml_params=None):

        if xml_node is not None and xml_params is not None:
            Reaction.__init__xml(self, xml_node, xml_params)

        else:

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

    def eval_score(self, score_dict):
        """
        Evalutes a reaction score.
        Returns 'nan' if there is no gene association
        """

        if self.gene_associations is None:
            return float('nan')
        else:
            return self.gene_associations.eval_score(score_dict)


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

    def eval_score(self, gene_scores):
        """
        Matches genes with scores
        Combines 'or' associations with 'max'
        Combines 'and' associations with 'min'
        """

        if self.type == 'and':
            return min([x.eval_score(gene_scores) for x in self.children])

        elif self.type == 'or':
            return max([x.eval_score(gene_scores) for x in self.children])

        elif self.type == 'gene':

            try:
                return gene_scores[self.gene.name]
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
