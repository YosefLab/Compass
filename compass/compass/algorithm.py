"""
Run the procedure for COMPASS
"""
from __future__ import print_function, division, absolute_import
import numpy as np
import pandas as pd
from tqdm import tqdm
from random import shuffle
import logging
import os
import sys

from .. import utils
from .. import models
from . import cache
from ..globals import BETA, EXCHANGE_LIMIT
import compass.global_state as global_state

import cplex

logger = logging.getLogger("compass")

__all__ = ['singleSampleCompass']


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
            - lambda, num_neighbors, symmetric_kernel, species,
              and_function, test_mode
    """

    if not os.path.isdir(directory):
        os.makedirs(directory)

    if os.path.exists(os.path.join(directory, 'success_token')):
        logger.info('success_token detected, results already calculated.')
        logger.info('COMPASS Completed Successfully')
        return

    model = models.init_model(model, species=args['species'],
                              exchange_limit=EXCHANGE_LIMIT,
                              media=media)

    logger.info("Running COMPASS on model: %s", model.name)

    if args['generate_cache']:
        cache.clear(model)

    # Build model into cplex problem
    problem = initialize_cplex_problem(model, args['num_threads'])

    # Only read this to get the number of samples and the sample name
    # Use nrows=1 so this is fast
    expression = pd.read_csv(data, sep='\t', index_col=0, nrows=1)
    sample_name = expression.columns[sample_index]

    logger.info("Processing Sample %i/%i: %s", sample_index,
                len(expression.columns), sample_name)
    global_state.set_current_cell_name(sample_name)
    # Run core compass algorithm

    # Evaluate reaction penalties
    logger.info("Evaluating Reaction Penalties...")
    reaction_penalties = pd.read_csv(
        args['penalties_file'], sep="\t", header=0,
        usecols=["Reaction", sample_name])

    reaction_penalties = reaction_penalties.set_index("Reaction").iloc[:, 0]

    if not args['no_reactions']:
        logger.info("Evaluating Reaction Scores...")
        reaction_scores = compass_reactions(
            model, problem, reaction_penalties,
            select_reactions=args['select_reactions'],
            TEST_MODE=args['test_mode'])

    #if user wants to calc reaction scores, but doesn't want to calc metabolite scores, calc only the exchange reactions
    logger.info("Evaluating Exchange/Secretion/Uptake Scores...")
    uptake_scores, secretion_scores, exchange_rxns = compass_exchange(
        model, problem, reaction_penalties,
        only_exchange=(not args['no_reactions']) and not args['calc_metabolites'],
        select_reactions=args['select_reactions'],
        TEST_MODE=args['test_mode'])

    # Copy valid uptake/secretion reaction fluxes from uptake/secretion
    #   results into reaction results
    if (not args['no_reactions']) or args['calc_metabolites']:
        for r_id in exchange_rxns:
            assert r_id in model.reactions
            assert r_id not in reaction_scores
            reaction_scores[r_id] = exchange_rxns[r_id]

    # Output results to file
    logger.info("Writing output files...")
    if not args['no_reactions']:
        reaction_scores = pd.Series(reaction_scores, name=sample_name).sort_index()
        reaction_scores.to_csv(os.path.join(directory, 'reactions.txt'),
                               sep="\t", header=True)

    if args['calc_metabolites']:
        uptake_scores = pd.Series(uptake_scores, name=sample_name).sort_index()
        secretion_scores = pd.Series(secretion_scores, name=sample_name).sort_index()

        uptake_scores.to_csv(os.path.join(directory, 'uptake.txt'),
                             sep="\t", header=True)
        secretion_scores.to_csv(os.path.join(directory, 'secretions.txt'),
                                sep="\t", header=True)

    if args['generate_cache']:
        logger.info(
            'Saving cache file for Model: {}, Media: {}'.format(
                model.name, model.media)
        )
        cache.save(model)

    # write success token
    with open(os.path.join(directory, 'success_token'), 'w') as fout:
        fout.write('Success!')

    logger.info('COMPASS Completed Successfully')


def compass_exchange(model, problem, reaction_penalties, select_reactions=None, only_exchange=False, TEST_MODE=False):
    """
    Iterates through metabolites, finding each's max
    uptake and secretion potentials. If only_exchange=True, does so only for exchange reactions.

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


    #populate the list of selected_reaction_ids - do this once outside of the loop
    if select_reactions:
        #assume this is a filename with one reaction per row, ignores unrecognized reactions
        if not os.path.exists(select_reactions):
            raise Exception("cannot find selected reactions subset file %s" % select_reactions)

        with open(select_reactions) as f:
            selected_reaction_ids = [line.strip() for line in f]


    all_names = set(problem.linear_constraints.get_names())

    for metabolite in tqdm(metabolites, file=sys.stderr):

        met_id = metabolite.id

        if met_id not in all_names:
            # This can happen if the metabolite does not participate
            # in any reaction. As a result, it won't be in any
            # constraints - happens in RECON2

            uptake_scores[met_id] = 0.0
            secretion_scores[met_id] = 0.0
            continue  # In test mode this always continues!

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

        #If user wants only exchange reaction - limit the reactions space through which we iterate
        if only_exchange:
            reactions = [x for x in reactions if x.is_exchange]

        if select_reactions:
            #r.id is a unidirectional identifier (ending with _pos or _neg suffix --> we remove it and compare to the undirected reaction id)
            reactions = [r for r in reactions if ((r.id)[:-4] in selected_reaction_ids)]


        for reaction in reactions:  # This only effectively loops over exchange reactions in fact.
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

        #if the selected_rxns or only_exchange options are used --> then we don't want to add reactions unless one of the pair already exists
        if(only_exchange or select_reactions) and (uptake_rxn is None) and (secretion_rxn is None):
            continue

        if (secretion_rxn is None):
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


        #if only exchange flag is set - don't add uptakes that do not exist
        if (uptake_rxn is None):
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

        global_state.set_current_reaction_id(secretion_rxn)
        value = solve_problem_wrapper(problem)
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

        global_state.set_current_reaction_id(uptake_rxn)
        value = solve_problem_wrapper(problem)
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


def compass_reactions(model, problem, reaction_penalties, select_reactions=None, TEST_MODE=False):
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

    if select_reactions:
        #assume this is a filename with one reaction per row, ignores unrecognized reactions
        if not os.path.exists(select_reactions):
            raise Exception("cannot find selected reactions subset file %s" % select_reactions)

        with open(select_reactions) as f:
            selected_reaction_ids = [line.strip() for line in f]

        #r.id is a unidirectional identifier (ending with _pos or _neg suffix --> we remove it and compare to the undirected reaction id)
        reactions = [r for r in reactions if ((r.id)[:-4] in selected_reaction_ids)]


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

            global_state.set_current_reaction_id(reaction.id)
            value = solve_problem_wrapper(problem)
            reaction_scores[reaction.id] = value

            # Remove Constraint
            problem.linear_constraints.delete('REACTION_OPT')

        # Restore limit of partner reaction to old state
        if partner_reaction is not None:
            partner_id = partner_reaction.id
            problem.variables.set_upper_bounds(partner_id, old_partner_ub)

    return reaction_scores


def initialize_cplex_problem(model, num_threads=1):
    # type: (compass.models.MetabolicModel)
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
    problem.parameters.threads.set(num_threads)

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

    global_state.set_current_reaction_id(rxn)
    rxn_max = solve_problem_wrapper(problem)

    # Save the result
    model_cache = cache.load(model)
    model_cache[rxn] = rxn_max

    return rxn_max


def solve_problem_wrapper(problem) -> float:
    r"""
    Only solve the problem if the reaction is selected for the cell. Else,
    skip the computation and return np.nan.
    """
    if global_state.current_reaction_is_selected_for_current_cell():
        problem.solve()
        return problem.solution.get_objective_value()
    else:
        return np.nan
