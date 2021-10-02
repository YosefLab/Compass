"""
For working with metabolic models
"""
from __future__ import print_function, division, absolute_import
import os
import json
import pandas as pd
from ..globals import MODEL_DIR
from math import isnan

# ----------------------------------------
# Functions for aggregation (and/or)
# ----------------------------------------


def min_w_nan(vals):
    """
    Min which propagates nan.
    Normal python 'min' ignores nan.
    """
    if any([isnan(x) for x in vals]):
        return float('nan')
    else:
        return min(vals)


def mean_nan_zero(vals):
    """
    Mean which treats 'nan' values as zeros
    """
    vals = [0 if isnan(x) else x for x in vals]
    return sum(vals) / len(vals)


def median_nan_zero(vals):
    """
    Median which treats 'nan' values as zeros
    """
    vals = [0 if isnan(x) else x for x in vals]
    vals = sorted(vals)
    if len(vals) % 2 == 1:
        middle_i = int((len(vals)-1)/2)
        return vals[middle_i]
    else:
        right_i = int(len(vals)/2)
        left_i = right_i-1
        return (vals[left_i] + vals[right_i])/2


def sum_wo_nan(vals):
    """
    Max which ignores nan.
    Normal python sum propagates nan.
    """
    vals = [x for x in vals if not isnan(x)]
    return sum(vals)


# ----------------------------------------
# Model class and related classes
# ----------------------------------------


class MetabolicModel(object):

    def __init__(self, name):
        self.name = name
        self.reactions = {}
        self.species = {}
        self.compartments = {}
        self.objectives = {}
        self._maximum_flux = None
        self.media = 'NoMedia'

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

    def getReactionExpression(self, expression,
                              and_function='min',
                              or_function=sum_wo_nan):
        # type: (pandas.Series) -> dict
        """
        Evaluates a score for every reaction, using the expression data.
        This is used for constraint/penalty generation.
        If a score cannot be computed, NaN is used to fill

        Returns a dict:
            key: reaction id (str)
            val: reaction score (float)
        """

        # resolve the AND function
        if and_function == 'min':
            and_function = min_w_nan
        elif and_function == 'mean':
            and_function = mean_nan_zero
        elif and_function == 'median':
            and_function = median_nan_zero
        else:
            raise ValueError("Invalid value for and_function: " +
                             str(and_function))

        score_dict = {}

        for r_id, reaction in self.reactions.items():
            score = reaction.eval_expression(
                expression, and_function, or_function)

            score_dict[r_id] = score

        return score_dict

    def limitExchangeReactions(self, limit):
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
                # Metabolites created as products - limit forward flux
                if reaction.upper_bound > limit:
                    reaction.upper_bound = limit

            elif len(reaction.products) == 0 and len(reaction.reactants) > 0:
                # Metabolites created as reactants - limit reverse flux
                if reaction.lower_bound < -1*limit:
                    reaction.lower_bound = -1*limit

            else:
                raise Exception("Should not occur")

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

            if reaction.is_pos_unidirectional:  # just ensure suffix and continue

                if (not reaction.id.endswith('_pos')) and \
                        (not reaction.id.endswith('_neg')):
                    reaction.id = reaction.id + "_pos"

                uni_reactions[reaction.id] = reaction

                continue

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

            # Only add reactions if they can carry flux
            if pos_reaction.upper_bound > 0:
                neg_reaction.reverse_reaction = pos_reaction
                uni_reactions[pos_reaction.id] = pos_reaction

            if neg_reaction.upper_bound > 0:
                pos_reaction.reverse_reaction = neg_reaction
                uni_reactions[neg_reaction.id] = neg_reaction

        self.reactions = uni_reactions

    def load_media(self, media_name):
        """
        Loads information in the media file and uses it to
        modify exchange reaction bounds

        Media files are stored in the model's directory under
        the `media` folder with names `<media_name>.json`

        Media file contains a JSON dict with keys corresponding to
        reaction IDs and values corresponding to reaction
        upper-bounds
        """

        media_file = media_name + '.json'
        media_file = os.path.join(MODEL_DIR, self.name, 'media', media_file)

        with open(media_file) as fin:
            media = json.load(fin)

        for rid, ub in media.items():
            self.reactions[rid].upper_bound = ub

        self.media = media_name

    def isoform_summing(self):
        """
        Removes instances where two isoforms of the same gene are summed/OR'd together
        """
        for reaction in self.reactions.values():
            reaction.gene_associations.isoform_summing()

    def _calc_max_flux(self):
        """
        Determines the max (absolute) flux of the model
        """
        max_flux = 0
        for reaction in self.reactions.values():
            max_flux = max(abs(reaction.lower_bound),
                           abs(reaction.upper_bound),
                           max_flux)

        self._maximum_flux = max_flux

    @property
    def maximum_flux(self):
        if self._maximum_flux is None:
            self._calc_max_flux()

        return self._maximum_flux

    @property
    def reaction_meta(self):
        out = pd.DataFrame()
        rows = [pd.Series(x.meta, name=x.id) for x in self.reactions.values()]
        out = out.append(rows)
        return out

    @property
    def species_meta(self):
        out = pd.DataFrame()
        rows = [pd.Series(x.meta, name=x.id) for x in self.species.values()]
        out = out.append(rows)
        return out

    def to_JSON(self):
        return json.dumps({
            'name': self.name,
            'reactions': self.reactions,
            'species': self.species,
            'compartments': self.compartments,
            'media': self.media
        }, indent=4, default=_json_default)


class Reaction(object):
    # Bounds
    # Reactants (and coefficients)
    # Products (and coefficients)
    # Gene Associations

    def __init__(self, from_reaction=None):

        if from_reaction is not None:  # Copy constructor

            self.id = from_reaction.id
            self.name = from_reaction.name
            self.subsystem = from_reaction.subsystem
            self.upper_bound = from_reaction.upper_bound
            self.lower_bound = from_reaction.lower_bound
            self.reactants = from_reaction.reactants.copy()
            self.products = from_reaction.products.copy()
            self.gene_associations = from_reaction.gene_associations
            self.reverse_reaction = from_reaction.reverse_reaction
            self.meta = from_reaction.meta

        else:  # Placeholders

            self.id = ""
            self.name = ""
            self.subsystem = ""
            self.upper_bound = float('nan')
            self.lower_bound = float('nan')
            self.reactants = {}
            self.products = {}
            self.gene_associations = None
            self.reverse_reaction = None
            self.meta = {}

    def eval_expression(self, expression, and_function, or_function):
        """
        Evalutes a reaction expression.
        Returns 'nan' if there is no gene association
        """

        if self.gene_associations is None:
            return float('nan')
        else:
            return self.gene_associations.eval_expression(expression, and_function, or_function)

    def list_genes(self):
        """
        Return all the genes associated with the reaction
        """
        if self.gene_associations is None:
            return []
        else:
            return list(self.gene_associations.list_genes())

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

    @property
    def is_pos_unidirectional(self):

        if (self.upper_bound > 0 and
                self.lower_bound == 0):
            return True
        else:
            return False

    def to_serializable(self):
        out = {
            'id': self.id,
            'name': self.name,
            'subsystem': self.subsystem,
            'upper_bound': self.upper_bound,
            'lower_bound': self.lower_bound,
            'reactants': self.reactants,
            'products': self.products,
            'meta': self.meta,
        }

        if self.reverse_reaction is not None:
            out['reverse_reaction'] = self.reverse_reaction.id
        else:
            out['reverse_reaction'] = ''

        # Gather a list of genes and their alt symbols
        def _get_genes(assoc):
            if assoc.type == 'gene':
                return [assoc.gene]
            else:
                result = []
                for ac in assoc.children:
                    result.extend(_get_genes(ac))
                return result

        gene_dict = {}
        if self.gene_associations is not None:
            genes = _get_genes(self.gene_associations)
            for gene in genes:
                gene_dict[gene.id] = gene

        out['genes'] = gene_dict

        out['gene_associations'] = self.gene_associations

        return out


class Species(object):

    def __init__(self):
            self.id = ""
            self.name = ""
            self.compartment = ""
            self.formula = ""
            self.meta = {}

    def to_serializable(self):
        return {
            'id': self.id,
            'name': self.name,
            'compartment': self.compartment,
            'formula': self.formula,
            'meta': self.meta
        }


class Compartment(object):

    def __init__(self, xml_node=None):
            self.id = ""
            self.name = ""

    def to_serializable(self):
        return {
            'id': self.id,
            'name': self.name
        }


class Association(object):

    def __init__(self):
            self.type = ''
            self.gene = None
            self.children = []

    def eval_expression(self, expression, and_function, or_function):
        """
        Matches genes with scores
        Combines 'or' associations with 'sum'
        Combines 'and' associations with 'min'
        """

        if self.type == 'and':
            return and_function(
                    [x.eval_expression(expression, and_function, or_function)
                     for x in self.children]
                    )

        elif self.type == 'or':
            return or_function(
                    [x.eval_expression(expression, and_function, or_function)
                     for x in self.children]
                    )

        elif self.type == 'gene':

            try:
                return self.gene.eval_expression(expression)
            except KeyError:
                return float('nan')

        else:
            raise Exception("Unknown Association Type")

    def list_genes(self):
        """
        Returns all the genes in the association
        """
        if self.type == 'gene':
            return {self.gene.name} | set(self.gene.alt_symbols)

        else:
            genes = set()

            for child in self.children:
                genes |= child.list_genes()

            return genes

    def isoform_summing(self):
        """
        Removes instances where two isoforms of the same gene are summed/OR'd together
        """
        if self.type == 'gene':
            return 

        elif self.type == 'or' or self.type == 'and':
            seen = set()
            children = []
            for child in self.children:
                if child.type == 'gene':
                    if child.gene.name not in seen:
                        seen.add(child.gene.name)
                        children.append(child)
                else:
                    child.isoform_summing()
                    children.append(child)
            self.children = children

    def __str__(self):
        return _print_node(self)

    def to_serializable(self):
        if self.type == 'gene':
            return {
                'type': self.type,
                'name': self.gene.name
            }
        else:
            return {
                'type': self.type,
                'children': [x.to_serializable() for x in self.children]
            }

class Gene(object):

    def __init__(self):
            self.id = ''
            self.name = ''
            self.alt_symbols = []

    def eval_expression(self, expression):
        """
        Get the expression for this gene in the expression vector

        Prefer match by name.  Alternately, match by alt_symbols
        """

        if self.name == '':
            return float('nan')

        if self.name in expression.index:
            return expression[self.name]

        # Average expression across found alt_symbols
        found_symbols = 0
        agg_expression = 0
        for symbol in self.alt_symbols:
            if symbol in expression.index:
                agg_expression += expression[symbol]
                found_symbols += 1

        if found_symbols > 0:
            return agg_expression / found_symbols

        return float('nan')

    def to_serializable(self):
        return {
            'id': self.id,
            'name': self.name,
            'alt_symbols': self.alt_symbols
        }


def _print_node(node, expression=None, indent=0):
    lines = []
    if expression is None:
        lines.append(" "*indent + node.type)
    else:
        lines.append(" ".join([
            " "*indent + node.type,
            node.eval_expression(expression, min_w_nan, sum_wo_nan)
        ]))

    if len(node.children) > 0:
        for x in node.children:
            lines.append(_print_node(x, expression, indent+4))
    if node.type == 'gene':
        lines.append(" ".join([
            " "*indent, node.gene.id,
            node.gene.name, str(node.gene.alt_symbols)]))

    return "\n".join(lines)


def _json_default(o):

    try:
        return o.to_serializable()
    except AttributeError:
        return o
