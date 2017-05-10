"""
Run the procedure for COMPASS
"""
from __future__ import print_function, division, absolute_import
import os
import json
import io
import pandas as pd
from tqdm import tqdm
from random import shuffle
import logging

logger = logging.getLogger("mflux")

from . import utils

import cplex

__all__ = ['run_compass']

from ..globals import RESOURCE_DIR
PREPROCESS_CACHE_DIR = os.path.join(RESOURCE_DIR, 'COMPASS')
if not os.path.isdir(PREPROCESS_CACHE_DIR):
    os.mkdir(PREPROCESS_CACHE_DIR)

BETA = 0.95  # Used to constrain model near optimal point
EXCHANGE_LIMIT = 1000.0  # Limit for exchange reactions


def run_compass(model, expression):
    # type: (mflux.models.MetabolicModel, pandas.DataFrame)
    """
    Runs COMPASS on many samples
    """

    if isinstance(expression, pd.Series):
        expression = pd.DataFrame(expression)

    logger.info("Running COMPASS on model: %s", model.name)

    # Limit exchange reactions
    model.limitExchangeReactions(limit=EXCHANGE_LIMIT)

    # Split fluxes into _pos / _neg
    model.make_unidirectional()

    # Build model into cplex problem
    problem = initialize_cplex_problem(model)

    # For each cell/sample
    # Eval expression scores
    # Eval reaction penalties
    # Eval Metabolite secretion/uptake penalties

    all_reaction_scores = {}
    all_uptake_scores = {}
    all_secretion_scores = {}

    for i, sample in enumerate(expression.columns):
        logger.info("Processing Sample %i/%i: %s", i,
                     len(expression.columns), sample)

        expression_data = expression[sample]

        reaction_penalties = eval_reaction_penalties(model, expression_data)

        reaction_scores = compass_reactions(
            model, problem, reaction_penalties)

        uptake_scores, secretion_scores = compass_exchange(
            model, problem, reaction_penalties)

        all_reaction_scores[sample] = pd.Series(reaction_scores)
        all_uptake_scores[sample] = pd.Series(uptake_scores)
        all_secretion_scores[sample] = pd.Series(secretion_scores)

    reaction_table = pd.DataFrame(all_reaction_scores)
    uptake_table = pd.DataFrame(all_uptake_scores)
    secretion_table = pd.DataFrame(all_secretion_scores)

    return reaction_table, uptake_table, secretion_table


def eval_reaction_penalties(model, expression_data):
    # type: (mflux.models.MetabolicModel, pandas.Series)
    """
    Determines reaction penalties, on the given model, for
    the given expression data
    """

    reaction_expression = model.getReactionExpression(expression_data)
    reaction_expression = pd.Series(reaction_expression)
    reaction_expression[reaction_expression < 0] = 0

    reaction_penalties = 1 / (1 + reaction_expression)

    reaction_penalties = reaction_penalties.dropna()

    return reaction_penalties


def compass_exchange(model, problem, reaction_penalties):
    """
    Iterates through metabolites, finding each's max
    uptake and secretion potentials.

    Holds each near its max uptake/secretion while minimizing
    penalty

    Returns the optimal penalty for uptake and secretion

    Returns
    -------
    uptake_scores: dict
        key: species_id
        value: minimum penalty achieved

    secretion_scores: dict
        key: species_id
        value: minimum penalty achieved
    """

    secretion_scores = {}
    uptake_scores = {}
    metabolites = list(model.species.values())
    shuffle(metabolites)


    all_names = set(problem.linear_constraints.get_names())

    for metabolite in tqdm(metabolites):

        met_id = metabolite.id

        if met_id not in all_names:
            # This can happen if the metabolite does not participate in any reaction
            # As a result, it won't be in any constraints - happened in RECON2

            uptake_scores[met_id] = 0.0
            secretion_scores[met_id] = 0.0
            continue

        # Rectify exchange reactions
        # Either find existing pos and neg exchange reactions
        # Or create new ones

        uptake_rxn = None
        extra_uptake_rxns = []
        secretion_rxn = None
        extra_secretion_rxns = []

        added_uptake = False     # Did we add an uptake reaction?
        added_secretion = False  #  "  "   "  "  secretion reaction? 

        # Metabolites represented by a constraint: get associated reactions
        sp = problem.linear_constraints.get_rows(met_id)
        rxn_ids = problem.variables.get_names(sp.ind)
        reactions = [model.reactions[x] for x in rxn_ids]

        for reaction in reactions:
            if reaction.is_exchange and met_id in reaction.products:
                if uptake_rxn is None:
                    uptake_rxn = reaction.id
                else:
                    extra_uptake_rxns.append(reaction.id)

            elif reaction.is_exchange and met_id in reaction.reactants:
                if secretion_rxn is None:
                    secretion_rxn = reaction.id
                else:
                    extra_secretion_rxns.append(reaction.id)

        if secretion_rxn is None:
            added_secretion = True
            secretion_rxn = met_id + "_SECRETION"

            # Add secretion reaction to the problem as a variable
            problem.variables.add(
                names=[secretion_rxn],
                ub=[model.maximum_flux],
                lb=[0.0],)

            # Add it to the metabolites constraint
            rxn_index = problem.variables.get_indices(secretion_rxn)
            sp.ind.append(rxn_index)
            sp.val.append(-1.0)

        if uptake_rxn is None:
            added_uptake = True
            uptake_rxn = met_id + "_UPTAKE"

            # Add uptake reaction to the problem as a variable
            problem.variables.add(
                names=[uptake_rxn],
                ub=[EXCHANGE_LIMIT],
                lb=[0.0],)

            # Add it to the metabolite's constraint
            rxn_index = problem.variables.get_indices(uptake_rxn)
            sp.ind.append(rxn_index)
            sp.val.append(1.0)

        # Modify the constraint in the problem
        #   e.g. Add the metabolites connections
        problem.linear_constraints.set_linear_components(met_id, sp)

        all_uptake = [uptake_rxn] + extra_uptake_rxns
        all_secretion = [secretion_rxn] + extra_secretion_rxns

        # -----------------
        # Optimal Secretion
        # -----------------

        # Close all uptake, storing their upper-bounds to restore later
        old_uptake_upper = {}
        for rxn_id in all_uptake:
            old_ub = problem.variables.get_upper_bounds(rxn_id)
            old_uptake_upper[rxn_id] = old_ub
            problem.variables.set_upper_bounds(rxn_id, 0.0)

        # Close extra secretion, storing upper-bounds to restore later
        old_secretion_upper = {}
        for rxn_id in extra_secretion_rxns:
            old_ub = problem.variables.get_upper_bounds(rxn_id)
            old_secretion_upper[rxn_id] = old_ub
            problem.variables.set_upper_bounds(rxn_id, 0.0)


        # Get max of secretion reaction
        secretion_max = maximize_reaction(model, problem, secretion_rxn)

        # Set contraint of max secretion to BETA*max
        problem.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=[secretion_rxn], val=[1.0])],
            senses=['R'],
            rhs=[BETA * secretion_max],
            names=['SECRETION_OPT'])

        # Find minimimum penalty
        utils.reset_objective(problem)
        problem.objective.set_linear(
            list(reaction_penalties.iteritems())
        )

        problem.objective.set_sense(problem.objective.sense.minimize)
        problem.solve()
        value = problem.solution.get_objective_value()
        secretion_scores[met_id] = value

        # Clear Secretion constraint
        problem.linear_constraints.delete('SECRETION_OPT')

        # Restore all uptake
        for rxn_id, old_ub in old_uptake_upper.items():
            problem.variables.set_upper_bounds(rxn_id, old_ub)

        # Restore extra secretion
        for rxn_id, old_ub in old_secretion_upper.items():
            problem.variables.set_upper_bounds(rxn_id, old_ub)

        # -----------------
        # Optimal Uptake
        # -----------------

        # Close extra uptake
        old_uptake_upper = {}
        for rxn_id in extra_uptake_rxns:
            old_ub = problem.variables.get_upper_bounds(rxn_id)
            old_uptake_upper[rxn_id] = old_ub
            problem.variables.set_upper_bounds(0.0)

        # Close all secretion
        old_secretion_upper = {}
        for rxn_id in all_secretion:
            old_ub = problem.variables.get_upper_bounds(rxn_id)
            old_secretion_upper[rxn_id] = old_ub
            problem.variables.set_upper_bounds(rxn_id, 0.0)

        # Get max of uptake reaction
        uptake_max = maximize_reaction(model, problem, uptake_rxn)

        # Set contraint of max uptake with BETA*max
        problem.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=[uptake_rxn], val=[1.0])],
            senses=['R'],
            rhs=[BETA * uptake_max],
            names=['UPTAKE_OPT'])

        # Find minimimum penalty
        utils.reset_objective(problem)
        problem.objective.set_linear(
            list(reaction_penalties.iteritems())
        )

        problem.objective.set_sense(problem.objective.sense.minimize)
        problem.solve()
        value = problem.solution.get_objective_value()
        uptake_scores[met_id] = value

        # Clear Secretion constraint
        problem.linear_constraints.delete('UPTAKE_OPT')

        # Restore extra uptake
        for rxn_id, old_ub in old_uptake_upper.items():
            problem.variables.set_upper_bounds(rxn_id, old_ub)

        # Restore all secretion
        for rxn_id, old_ub in old_secretion_upper.items():
            problem.variables.set_upper_bounds(rxn_id, old_ub)

        # Remove added uptake and secretion reactions
        if added_uptake:
            problem.variables.delete(uptake_rxn)

        if added_secretion:
            problem.variables.delete(secretion_rxn)

    return uptake_scores, secretion_scores


def compass_reactions(model, problem, reaction_penalties):
    """
    Iterates through reactions, holding each near
    its max value while minimizing penalty.

    Minimum overall penalty returned for each reaction

    Returns
    -------
    reaction_scores: dict
        key: reaction id
        value: minimum penalty achieved
    """

    # Create overall penalty as the objective function to minimize
    utils.reset_objective(problem)
    problem.objective.set_linear(
        list(reaction_penalties.iteritems())
    )

    # Iterate through Reactions
    reaction_scores = {}

    reactions = list(model.reactions.values())
    shuffle(reactions)
    for reaction in reactions:

        if reaction.is_exchange:
            continue

        partner_reaction = reaction.reverse_reaction

        # Set partner reaction upper-limit to 0 in problem
        # Store old limit for later to restore
        if partner_reaction is not None:
            partner_id = partner_reaction.id
            old_partner_ub = problem.variables.get_upper_bounds(partner_id)
            problem.variables.set_upper_bounds(partner_id, 0.0)

        r_max = maximize_reaction(model, problem, reaction.id)

        # If Reaction can't carry flux anyways, just continue
        if r_max == 0:
            reaction_scores[reaction.id] = 0

        else:

            problem.linear_constraints.add(
                lin_expr=[cplex.SparsePair(ind=[reaction.id], val=[1.0])],
                senses=['R'],
                rhs=[BETA * r_max],
                names=['REACTION_OPT'])

            # Minimize Penalty
            utils.reset_objective(problem)
            problem.objective.set_linear(
                list(reaction_penalties.iteritems())
            )

            problem.objective.set_sense(problem.objective.sense.minimize)
            problem.solve()
            value = problem.solution.get_objective_value()
            reaction_scores[reaction.id] = value

            # Remove Constraint
            problem.linear_constraints.delete('REACTION_OPT')

        # Restore limit of partner reaction to old state
        if partner_reaction is not None:
            partner_id = partner_reaction.id
            problem.variables.set_upper_bounds(partner_id, old_partner_ub)

    return reaction_scores


def initialize_cplex_problem(model):
    # type: (mflux.models.MetabolicModel)
    """
    Builds and returns a cplex problem representing our metabolic model

    Limits exchange reactions and makes all reactions unidirectional
    by splitting into two components
    """

    # Create the Problem first
    # Easier to modify existing problem and re-solve
    problem = cplex.Cplex()
    problem.set_log_stream(None)  # Suppress output
    problem.set_error_stream(None)  # Suppress errors
    problem.set_warning_stream(None)  # Suppress Warnings
    problem.set_results_stream(None)  # Suppress results to output

    # Set Parameters for the Cplex solver
    problem.parameters.emphasis.numerical.set(True)

    # Add variables
    reactions = list(model.reactions.values())
    problem.variables.add(
        names=[x.id for x in reactions],
        ub=[x.upper_bound for x in reactions],
        lb=[x.lower_bound for x in reactions],)

    # Add constraints

    # Add stoichiometry constraints
    c_lin_expr, c_senses, c_rhs, c_names = (
        utils.get_steadystate_constraints(model))

    problem.linear_constraints.add(
        lin_expr=c_lin_expr,
        senses=c_senses,
        rhs=c_rhs,
        names=c_names)

    return problem


def maximize_reaction(model, problem, rxn, use_cache=True):
    """Maximizes the current reaction in the problem
    Attempts to retrieve the value from cache if its in cache
    """

    # Load from cache if it exists and return
    if use_cache:
        cache = _load_cache(model)
        if rxn in cache:
            return cache[rxn]

    # Maximize the reaction
    utils.reset_objective(problem)
    problem.objective.set_linear(rxn, 1.0)
    problem.objective.set_sense(problem.objective.sense.maximize)

    problem.solve()
    rxn_max = problem.solution.get_objective_value()

    # Save the result
    cache = _load_cache(model)
    cache[rxn] = rxn_max

    return rxn_max


_cache = {}
def _load_cache(model):
    global _cache

    if model.name not in _cache:

        cache_file = os.path.join(PREPROCESS_CACHE_DIR,
                                  model.name + ".preprocess")

        if os.path.exists(cache_file):
            with open(cache_file) as fin:
                out = json.load(fin)

            _cache[model.name] = out

        else:
            _cache[model.name] = {}

    return _cache[model.name]


def _save_cache(model):
    global _cache

    cache_data = _cache[model.name]

    cache_file = os.path.join(PREPROCESS_CACHE_DIR,
                              model.name + ".preprocess")

    with open(cache_file, 'w') as fout:
        json.dump(cache_data, fout, indent=1)


def _clear_cache(model):
    global _cache

    _cache[model.name] = {}

    _save_cache(model)
