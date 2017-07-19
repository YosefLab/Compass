"""
For working with metabolic models
"""
from __future__ import print_function, division
import os
import libsbml
import json
import re
import pandas as pd
from .globals import RESOURCE_DIR
from math import isnan
from six import string_types

MODEL_DIR = os.path.join(RESOURCE_DIR, 'Metabolic Models')

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
# Loading models from either XML or MATLAB outputs
# ----------------------------------------


def load_metabolic_model(model_name, species='homo_sapiens'):
    """
    Loads the metabolic model from `file_name`, returning a Model object
    """

    if model_name.endswith('_xml'):
        model = load_metabolic_model_xml(model_name)
    elif model_name.endswith('_mat'):
        model = load_metabolic_model_matlab(model_name, species)
    else:
        raise NotImplementedError(
            "Can only handle .xml files or _mat directories")

    return model


def load_metabolic_model_xml(model_name):

    model_dir = os.path.join(MODEL_DIR, model_name)
    model_file = [x for x in os.path.listdir(model_dir)
                  if x.lower().endswith('.xml')][0]
    full_path = os.path.join(model_dir, model_file)
    sbmlDocument = libsbml.readSBMLFromFile(full_path)
    xml_model = sbmlDocument.model

    # Load Model Parameters
    xml_params = {}
    for pp in xml_model.getListOfParameters():
        key = pp.getId()
        val = pp.getValue()
        xml_params[key] = val

    modelx = MetabolicModel(model_name)

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
        species = Species(xml_node=met)
        modelx.species[species.id] = species

    return modelx


def load_metabolic_model_matlab(model_name, species):
    """
    model_name: str
        Name of the folder containing the model
    species: str
        Species name.  either 'homo_sapiens' or 'mus_musculus'
    """

    # First load Genes
    top_dir = os.path.join(MODEL_DIR, model_name)
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
        rxnMeta = pd.read_table(fname, index_col=0, quoting=3)

        for i, reaction in enumerate(reactions):
            reaction.meta = rxnMeta.loc[reaction.id].to_dict()

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
        metMeta = pd.read_table(fname, index_col=0, quoting=3)

        for i, met in enumerate(species):
            met.meta = metMeta.loc[met.id].to_dict()

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

    model = MetabolicModel(name)
    model.reactions = reactions
    model.species = species

    return model

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
        self.media = 'None'

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

    def __init__xml(self, xml_node, xml_params):
        """
        Build the reaction from the xml node

        Needs xml_params for bounds
        """
        self.id = xml_node.getId()
        self.name = xml_node.getName()
        self.subsystem = ""
        self.reverse_reaction = None
        self.meta = {}

        # Lower and upper bounds

        fbcrr = xml_node.getPlugin('fbc')
        ub = fbcrr.getUpperFluxBound()
        if isinstance(ub, string_types):
            ub = xml_params[ub]

        lb = fbcrr.getLowerFluxBound()
        if isinstance(lb, string_types):
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


class Species(object):

    def __init__(self, xml_node=None):
        if xml_node is not None:
            Species.__init__xml(self, xml_node)
        else:
            self.id = ""
            self.name = ""
            self.compartment = ""
            self.formula = ""
            self.meta = {}

    def __init__xml(self, xml_node):

        self.id = xml_node.getId()
        self.name = xml_node.getName()
        self.compartment = xml_node.getCompartment()
        self.formula = xml_node.getPlugin('fbc') \
            .getChemicalFormula()
        self.meta = {}


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


class Gene(object):

    def __init__(self, xml_node=None):
        if xml_node is not None:
            Gene.__init__xml(self, xml_node)
        else:
            self.id = ''
            self.name = ''
            self.alt_symbols = []

    def __init__xml(self, xml_node):

        self.id = xml_node.getId()
        self.name = xml_node.getName().upper()
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

# ----------------------------------------
# Utility functions
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



def _print_node(node, expression=None, indent=0):
    if expression is None:
        print(" "*indent + node.type)
    else:
        print(" "*indent + node.type, node.eval_expression(expression, min_w_nan, sum_wo_nan))

    if len(node.children) > 0:
        for x in node.children:
            _print_node(x, expression, indent+4)
    if node.type == 'gene':
        print(" "*indent, node.gene.id, node.gene.name, node.gene.alt_symbols)
