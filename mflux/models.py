"""
For working with metabolic models
"""
from __future__ import print_function, division
import os
import libsbml
import pandas as pd
from .globals import RESOURCE_DIR

MODEL_DIR = os.path.join(RESOURCE_DIR, 'Metabolic Models')


def load_metabolic_model(file_name):
    """
    Loads the metabolic model from `file_name`, returning a Model object
    """

    full_path = os.path.join(MODEL_DIR, file_name)
    sbmlDocument = libsbml.readSBMLFromFile(full_path)

    model = MetabolicModel(sbmlDocument)

    return model


class MetabolicModel(object):

    def __init__(self, sbmlDocument):
        self.sbmlDocument = sbmlDocument  # If you don't save this, get seg faults for some reason
        self.model = sbmlDocument.model
        self.fbmodel = self.model.getPlugin('fbc')

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

        reaction_list = self.model.getListOfReactions()
        score_dict = {}
        for reaction in reaction_list:
            fbc_reaction = reaction.getPlugin('fbc')

            gpa = fbc_reaction.getGeneProductAssociation()

            if gpa is None:
                score_dict[reaction.getId()] = float('nan')
            else:
                score = eval_Association(self.fbmodel, gpa.getAssociation(), expression)
                score_dict[reaction.getId()] = score

        return score_dict

    def getReactions(self):
        """
        Returns a list of reaction id's in the MetabolicModel
        """
        r_ids = []
        for rr in self.model.getListOfReactions():
            r_ids.append(rr.getId())

        return r_ids

    def getReactionBounds(self):
        """
        Returns two dicts, each mapping reaction id -> bound

        Returns:
            lower_bounds
            upper_bounds
        """
        # Load parameters
        params = {}
        for pp in self.model.getListOfParameters():
            key = pp.getId()
            val = pp.getValue()
            params[key] = val

        lbounds = {}
        ubounds = {}
        for rr in self.model.getListOfReactions():
            fbcrr = rr.getPlugin('fbc')
            ub = fbcrr.getUpperFluxBound()
            if isinstance(ub, str):
                ub = params[ub]

            lb = fbcrr.getLowerFluxBound()
            if isinstance(lb, str):
                lb = params[lb]

            reactionId = rr.getId()

            lbounds[reactionId] = lb
            ubounds[reactionId] = ub

        return lbounds, ubounds

    def getSMAT(self):
        """
        Returns a sparse form of the s-matrix

        result is a dict
            key: metabolite (species) id
            value: list of 2-tuples (reaction_id, coefficient)

        coefficient is positive if metabolite is produced in the reaction, negative if consumed
        """
        s_mat = {}
        for rr in self.model.getListOfReactions():

            reaction_id = rr.getId()

            # reactants
            for sr in rr.getListOfReactants():

                metabolite = sr.getSpecies()
                coefficient = sr.getStoichiometry()
                coefficient *= -1

                if metabolite not in s_mat:
                    s_mat[metabolite] = []

                s_mat[metabolite].append((reaction_id, coefficient))

            # products
            for sr in rr.getListOfProducts():

                metabolite = sr.getSpecies()
                coefficient = sr.getStoichiometry()

                if metabolite not in s_mat:
                    s_mat[metabolite] = []

                s_mat[metabolite].append((reaction_id, coefficient))

        return s_mat


def eval_Association(fbmodel, gpa, expression):
    # type: (libsbml.FbcModelPlugin, libsbml.FbcAssociation, pandas.Series) -> float
    """
    Parses the relationships in a FbcAssociation and evaluates the resulting expression value

    For an AND association, return the minimum
    For an OR association, return the max

    Currently, missing genes are treated as 0's.
    """

    if isinstance(gpa, libsbml.FbcOr):
        vals = [eval_Association(fbmodel, x, expression) for x in gpa.getListOfAssociations()]
        return max(vals)

    elif isinstance(gpa, libsbml.FbcAnd):
        vals = [eval_Association(fbmodel, x, expression) for x in gpa.getListOfAssociations()]
        return min(vals)

    elif isinstance(gpa, libsbml.GeneProductRef):
        gp_str = gpa.getGeneProduct()
        gp = fbmodel.getGeneProduct(gp_str)
        gene_name = gp.getName()

        try:
            return expression[gene_name]
        except KeyError:
            return 0.0

    else:
        raise ValueError("Unknown input type: " + type(gpa))
