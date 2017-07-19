"""
Run the procedure for COMPASS
"""
from __future__ import print_function, division, absolute_import
import pandas as pd
from tqdm import tqdm
from random import shuffle
import logging
import os
import sys

from .. import utils
from ..globals import NUM_THREADS, TEST_MODE
from .. import models
from . import cache
from . import penalties

import cplex

logger = logging.getLogger("mflux")

__all__ = ['singleSampleCompass']

BETA = 0.95  # Used to constrain model near optimal point
EXCHANGE_LIMIT = 1.0  # Limit for exchange reactions


def singleSampleCompass(data, model, media, directory, sample_index, args):
    """
    Run Compass on a single column of data

    Parameters
    ==========
    data : str
       Full path to data file

    model : str
        Name of metabolic model to use

    media : str or None
        Name of media to use

    directory : str
        Where to store results and log info.  Is created if it doesn't exist.

    sample_index : int
        Which sample to run on

    args : dict
        More keyword arguments
    """

    if not os.path.isdir(directory):
        os.makedirs(directory)

    # Unpack extra arguments
    lambda_ = args['lambda']
    perplexity = args['perplexity']
    symmetric_kernel = args['symmetric_kernel']

    expression = pd.read_table(data, index_col=0)
    expression.index = expression.index.str.upper()  # Gene names to upper
    model = models.load_metabolic_model(model, args['species'])
    logger.info("Running COMPASS on model: %s", model.name)

    # Limit exchange reactions
    model.limitExchangeReactions(limit=EXCHANGE_LIMIT)

    # Split fluxes into _pos / _neg
    model.make_unidirectional()

    if media is not None:
        model.load_media(media)

    # Build model into cplex problem
    problem = initialize_cplex_problem(model)

    expression_data = expression.iloc[:, sample_index]
    sample_name = expression.columns[sample_index]

    logger.info("Processing Sample %i/%i: %s", sample_index,
                len(expression.columns), sample_name)

    # Run core compass algorithm

    # Evaluate reaction penalties
    logger.info("Evaluating Reaction Penalties...")

    if lambda_ == 0:
        reaction_penalties = penalties.eval_reaction_penalties(
            model, expression_data,
            and_function=args['and_function'])
    else:
        reaction_penalties = penalties.eval_reaction_penalties_shared(
            model, expression, sample_index, lambda_, perplexity=perplexity,
            symmetric_kernel=symmetric_kernel,
            and_function=args['and_function'])

    logger.info("Evaluating Reaction Scores...")
    reaction_scores = compass_reactions(
        model, problem, reaction_penalties)

    logger.info("Evaluating Secretion/Uptake Scores...")
    uptake_scores, secretion_scores, exchange_rxns = compass_exchange(
        model, problem, reaction_penalties)

    # Copy valid uptake/secretion reaction fluxes from uptake/secretion
    #   results into reaction results

    for r_id in exchange_rxns:
        assert r_id in model.reactions
        assert r_id not in reaction_scores
        reaction_scores[r_id] = exchange_rxns[r_id]

    # Output results to file

    reaction_scores = pd.Series(reaction_scores, name=sample_name)
    uptake_scores = pd.Series(uptake_scores, name=sample_name)
    secretion_scores = pd.Series(secretion_scores, name=sample_name)

    reaction_scores.to_csv(os.path.join(directory, 'reactions.txt'),
                           sep="\t", header=True)
    uptake_scores.to_csv(os.path.join(directory, 'uptake.txt'),
                         sep="\t", header=True)
    secretion_scores.to_csv(os.path.join(directory, 'secretions.txt'),
                            sep="\t", header=True)

    logger.info('COMPASS Completed Successfully')


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

    exchange_rxns: dict
        Separate storage for exchange reactions.  These
        are skipped in the reaction loop.
        key: rxn_id
        value: minimum penalty achieved
    """

    secretion_scores = {}
    uptake_scores = {}
    exchange_rxns = {}
    metabolites = list(model.species.values())
    if TEST_MODE:
        metabolites = metabolites[0:50]
    shuffle(metabolites)

    all_names = set(problem.linear_constraints.get_names())

    for metabolite in tqdm(metabolites, file=sys.stderr):

        met_id = metabolite.id

        if met_id not in all_names:
            # This can happen if the metabolite does not participate
            # in any reaction. As a result, it won't be in any
            # constraints - happens in RECON2

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
        added_secretion = False  # "   "   "  "  secretion reaction?

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
        if(problem.objective.get_name() != 'reaction_penalties'):
            utils.reset_objective(problem)
            problem.objective.set_linear(
                list(reaction_penalties.iteritems())
            )
            problem.objective.set_name('reaction_penalties')
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
        if(problem.objective.get_name() != 'reaction_penalties'):
            utils.reset_objective(problem)
            problem.objective.set_linear(
                list(reaction_penalties.iteritems())
            )
            problem.objective.set_name('reaction_penalties')
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
        else:
            for rxn_id in all_uptake:
                exchange_rxns[rxn_id] = uptake_scores[met_id]

        if added_secretion:
            problem.variables.delete(secretion_rxn)
        else:
            for rxn_id in all_secretion:
                exchange_rxns[rxn_id] = secretion_scores[met_id]

    return uptake_scores, secretion_scores, exchange_rxns


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

    # Iterate through Reactions
    reaction_scores = {}

    reactions = list(model.reactions.values())
    if TEST_MODE:
        reactions = reactions[0:100]
    shuffle(reactions)

    for reaction in tqdm(reactions, file=sys.stderr):

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
            if(problem.objective.get_name() != 'reaction_penalties'):
                utils.reset_objective(problem)
                problem.objective.set_linear(
                    list(reaction_penalties.iteritems())
                )
                problem.objective.set_name('reaction_penalties')
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
    problem.parameters.threads.set(NUM_THREADS)

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

    # Initialize the objective
    utils.reset_objective(problem)

    return problem


def maximize_reaction(model, problem, rxn, use_cache=True):
    """Maximizes the current reaction in the problem
    Attempts to retrieve the value from cache if its in cache
    """

    # Load from cache if it exists and return
    if use_cache:
        model_cache = cache.load(model)
        if rxn in model_cache:
            return model_cache[rxn]

    # Maximize the reaction
    utils.reset_objective(problem)
    problem.objective.set_linear(rxn, 1.0)
    problem.objective.set_name(rxn)
    problem.objective.set_sense(problem.objective.sense.maximize)

    problem.solve()
    rxn_max = problem.solution.get_objective_value()

    # Save the result
    model_cache = cache.load(model)
    model_cache[rxn] = rxn_max

    return rxn_max
