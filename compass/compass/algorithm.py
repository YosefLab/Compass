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
import time
import timeit
import numpy as np
from .. import utils
from .. import models
from . import cache
from ..globals import BETA, EXCHANGE_LIMIT

import cplex

logger = logging.getLogger("compass")

__all__ = ['singleSampleCompass']


def singleSampleCompass(data, model, media, directory, sample_index, args):
    """
    Run Compass on a single column of data

    Parameters
    ==========
    data : list
       Full path to data file(s)

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
              and_function, test_mode, detailed_perf
    """
    if not os.path.isdir(directory) and directory != '/dev/null':
        os.makedirs(directory)

    if os.path.exists(os.path.join(directory, 'success_token')):
        logger.info('success_token detected, results already calculated.')
        logger.info('COMPASS Completed Successfully')
        return

    if args['save_argmaxes']:
        args['save_argmaxes_dir'] = directory
    else:
        args['save_argmaxes_dir'] = None

    model = models.init_model(model, species=args['species'],
                              exchange_limit=EXCHANGE_LIMIT,
                              media=media)

    if args['glucose']:
        model.reactions['EX_glc(e)_neg'].upper_bound = args['glucose']
        model.reactions['GLCt1r_pos'].upper_bound = args['glucose']

    logger.info("Running COMPASS on model: %s", model.name)

    perf_log = None
    if args['detailed_perf']:
        cols = ['order','max rxn time', 'max rxn method', 'cached', 'min penalty time', 
        'min penalty method', 'min penalty sensitvivity', 'kappa']
        perf_log = {c:{} for c in cols}

    if args['generate_cache']:
        cache.clear(model) #TBD add media specifier here too

    # Build model into cplex problem
    problem = initialize_cplex_problem(model, args['num_threads'], args['lpmethod'], args['advance'])

    # Only read this to get the number of samples and the sample name
    # Use nrows=1 so this is fast
    expression = utils.read_data(data)#pd.read_csv(data, sep='\t', index_col=0, nrows=1) #TBD add similar function for counting mtx stuff
    sample_name = str(expression.columns[sample_index])

    logger.info("Processing Sample %i/%i: %s", sample_index,
                len(expression.columns), sample_name)

    # Run core compass algorithm

    # Evaluate reaction penalties
    penalty_start = timeit.default_timer() #
    logger.info("Evaluating Reaction Penalties...")
    reaction_penalties = pd.read_csv(
        args['penalties_file'], sep="\t", header=0,
        usecols=["Reaction", sample_name])

    reaction_penalties = reaction_penalties.set_index("Reaction").iloc[:, 0]
    penalty_elapsed = timeit.default_timer() - penalty_start


    react_start = timeit.default_timer()
    if not args['no_reactions']:
        logger.info("Evaluating Reaction Scores...")
        reaction_scores = compass_reactions(
            model, problem, reaction_penalties,
            perf_log=perf_log, args=args)
    react_elapsed = timeit.default_timer() - react_start
    

    #if user wants to calc reaction scores, but doesn't want to calc metabolite scores, calc only the exchange reactions
    logger.info("Evaluating Exchange/Secretion/Uptake Scores...")
    exchange_start = timeit.default_timer()
    uptake_scores, secretion_scores, exchange_rxns = compass_exchange(
        model, problem, reaction_penalties,
        only_exchange=(not args['no_reactions']) and not args['calc_metabolites'],
        perf_log=perf_log, args=args)
    exchange_elapsed = timeit.default_timer() - exchange_start 

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

    if args['generate_cache'] or cache.is_new_cache(model):
        logger.info(
            'Saving cache file for Model: {}, Media: {}'.format(
                model.name, model.media)
        )
        cache.save(model)

    # write success token
    with open(os.path.join(directory, 'success_token'), 'w') as fout:
        fout.write('Success!')
    logger.info("Compass Penalty Time: ", penalty_elapsed)
    if not args['no_reactions']:
        logger.info("Compass Reaction Time:", react_elapsed,
        "Average time:", react_elapsed/len(reaction_scores))
        logger.info("Processed", len(reaction_scores), "reactions")
    logger.info("Compass Exchange Time:", exchange_elapsed, 
        "Average time:", exchange_elapsed/(len(secretion_scores)+len(uptake_scores)))
    logger.info("Processed ", len(uptake_scores), " uptake reactions")
    logger.info("Processed ", len(secretion_scores), " secretion reactions")

    if perf_log is not None:
        perf_log = pd.DataFrame(perf_log)
        perf_log.to_csv(os.path.join(directory, "compass_performance_log.csv"))
        logger.info("Saved detailed performance log")
    
    logger.info('COMPASS Completed Successfully')
    
def read_selected_reactions(select_reactions, select_subsystems, model):
    selected_reaction_ids = []
    if select_reactions:
        if not os.path.exists(select_reactions):
            raise Exception("cannot find selected reactions subset file %s" % select_reactions)
        with open(select_reactions) as f:
            selected_reaction_ids += [line.strip() for line in f]
    if select_subsystems:
        if not os.path.exists(select_subsystems):
            raise Exception("cannot find selected reactions subset file %s" % select_subsystems)
        with open(select_subsystems) as f:
            selected_subsystems_ids = [line.strip() for line in f]
        subsys = {}
        for rxn in model.reactions.values():
            if rxn.subsystem in selected_subsystems_ids:
                selected_reaction_ids += [rxn.id]
    return [str(s) for s in selected_reaction_ids]



def compass_exchange(model, problem, reaction_penalties, only_exchange=False, perf_log=None, args = None):
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
    if args['test_mode']:
        metabolites = metabolites[0:50]

    #populate the list of selected_reaction_ids - do this once outside of the loop
    if args['select_reactions'] or args['select_subsystems']:
        selected_reaction_ids = read_selected_reactions(args['select_reactions'], args['select_subsystems'], model)


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

        #If user wants only exchange reaction - limit the reactions space through which we iterate
        if only_exchange:
            reactions = [x for x in reactions if x.is_exchange]

        if args['select_reactions'] or args['select_subsystems']:
            #r.id is a unidirectional identifier (ending with _pos or _neg suffix --> we remove it and compare to the undirected reaction id)
            reactions = [r for r in reactions if ((r.id)[:-4] in selected_reaction_ids or str(r.id) in selected_reaction_ids)]


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

        #if the selected_rxns or only_exchange options are used --> then we don't want to add reactions unless one of the pair already exists
        if(only_exchange or args['select_reactions']) and (uptake_rxn is None) and (secretion_rxn is None):
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
        secretion_max = maximize_reaction(model, problem, secretion_rxn, use_cache=(args['glucose'] is None) , perf_log=perf_log)

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

        if perf_log is not None:
            start_time = timeit.default_timer()#time.perf_counter() #Not in python2.7
        problem.solve()
        if perf_log is not None:
            perf_log['min penalty time'][secretion_rxn] = timeit.default_timer() - start_time #time.perf_counter() - start_time #Not in python2.7
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
        uptake_max = maximize_reaction(model, problem, uptake_rxn, perf_log=perf_log)

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

        if perf_log is not None:
            start_time = timeit.default_timer() #time.perf_counter() #Not in python2.7
        problem.solve()
        if perf_log is not None:
            perf_log['min penalty time'][uptake_rxn] = timeit.default_timer() - start_time #time.perf_counter() - start_time #Not in python2.7
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


def compass_reactions(model, problem, reaction_penalties, perf_log=None, args = None):
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
    model_cache = cache.load(model)

    if args['test_mode']:
        reactions = reactions[0:100]

    if args['select_reactions'] or args['select_subsystems']:
        selected_reaction_ids = read_selected_reactions(args['select_reactions'], args['select_subsystems'], model)
        #r.id is a unidirectional identifier (ending with _pos or _neg suffix --> we remove it and compare to the undirected reaction id)
        reactions = [r for r in reactions if (str(r.id)[:-4] in selected_reaction_ids or str(r.id) in selected_reaction_ids)]

    
    if args['save_argmaxes']:
        argmaxes_order = []
        argmaxes = []

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

        
        r_max = maximize_reaction(model, problem, reaction.id, use_cache =(args['glucose'] is None), perf_log=perf_log)
        

        # If Reaction can't carry flux anyways, just continue
        if r_max == 0:
            reaction_scores[reaction.id] = 0
            if perf_log is not None:
               perf_log['min penalty time'][reaction.id] = 0
               #perf_log['blocked'][reaction.id] = True

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
            
            if perf_log is not None:
                #perf_log['blocked'][reaction.id] = False
                start_time = timeit.default_timer() #time.perf_counter() #Not in python2.7
            
            problem.solve()
            if perf_log is not None:
                perf_log['min penalty time'][reaction.id] = timeit.default_timer() - start_time #time.perf_counter() - start_time #Not in python2.7
                perf_log['min penalty method'][reaction.id] = problem.solution.get_method()
                perf_log['min penalty sensitvivity'][reaction.id] = problem.solution.sensitivity.objective(reaction.id)
                if hasattr(problem.solution.get_quality_metrics(),'kappa'):
                   perf_log['kappa'][reaction.id] = problem.solution.get_quality_metrics().kappa
            
            if args['save_argmaxes']:
                argmaxes.append(np.array(problem.solution.get_values()))
                argmaxes_order.append(reaction.id)

            value = problem.solution.get_objective_value()
            reaction_scores[reaction.id] = value

            # Remove Constraint
            problem.linear_constraints.delete('REACTION_OPT')

        # Restore limit of partner reaction to old state
        if partner_reaction is not None:
            partner_id = partner_reaction.id
            problem.variables.set_upper_bounds(partner_id, old_partner_ub)

    if args['save_argmaxes']:
        argmaxes = np.vstack(argmaxes)
        np.save(os.path.join(args['save_argmaxes_dir'],'argmaxes.npy'), argmaxes)
        argmaxes_order = np.array(argmaxes_order)
        np.save(os.path.join(args['save_argmaxes_dir'],'argmaxes_order.npy'), argmaxes_order)

    return reaction_scores


def initialize_cplex_problem(model, num_threads=1, lpmethod=0, adv=2):
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
    problem.parameters.preprocessing.reduce.set(3) #Turning on primal and dual preprocessing also enables some reoptimization features
    problem.parameters.advance.set(adv) #Will presolve advanced basis again
    problem.parameters.barrier.convergetol.set(1e-12) #default is 1e-8, minimum is 1e-12.
    problem.parameters.simplex.tolerances.optimality.set(1e-9) #default 1e-6, minimum is 1e-9
    problem.parameters.lpmethod.set(lpmethod) #default lets CPLEX choose the method

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


def maximize_reaction(model, problem, rxn, use_cache=True, perf_log=None):
    """Maximizes the current reaction in the problem
    Attempts to retrieve the value from cache if its in cache
    """
    
    if perf_log is not None:
        start_time = timeit.default_timer() #time.perf_counter() #Not in python2.7
        perf_log['order'][rxn] = len(perf_log['order'])

    # Load from cache if it exists and return
    if use_cache:
        model_cache = cache.load(model)
        if rxn in model_cache:
            if perf_log is not None:
                perf_log['cached'][rxn] = True
                perf_log['max rxn time'][rxn] = timeit.default_timer() - start_time #time.perf_counter() - start_time #Not in python2.7
            return model_cache[rxn]

    # Maximize the reaction
    utils.reset_objective(problem)
    problem.objective.set_linear(rxn, 1.0)
    problem.objective.set_name(str(rxn))
    problem.objective.set_sense(problem.objective.sense.maximize)

    problem.solve()
    rxn_max = problem.solution.get_objective_value()

    # Save the result
    model_cache = cache.load(model)
    model_cache[rxn] = rxn_max

    if perf_log is not None:
        perf_log['cached'][rxn] = False
        perf_log['max rxn time'][rxn] = timeit.default_timer() - start_time #time.perf_counter() - start_time #Not in python2.7
        perf_log['max rxn method'][rxn] = problem.solution.get_method()

    return rxn_max

def maximize_reaction_range(start_stop, args):
    """
    Maximizes a range of reactions from start_stop=(start, stop).
    args must be a dict with keys 'model', 'species', 'media'

    This is to reduce overhead compared to initializing the model and problem each time
    The model and problem cannot be passed to partial because they are not pickleable
    and pool requires the function to be pickleable.

    Returns
    -------
    sub_cache : dict
        key : reaction id
        value : maximum flux for reaction
    """
    #make a sub cache for each thread to write into
    sub_cache = {}
    model = models.init_model(args['model'], species=args['species'],
                              exchange_limit=EXCHANGE_LIMIT,
                              media=args['media'])
    problem = initialize_cplex_problem(model, args['num_threads'], args['lpmethod'])

    #sort by id to ensure consistency across threads
    reactions = sorted(list(model.reactions.values()), key=lambda r:r.id)[start_stop[0]:start_stop[1]]
    for reaction in tqdm(reactions, file=sys.stderr):
        #if reaction.is_exchange:
        #    continue
        partner_reaction = reaction.reverse_reaction
        # Set partner reaction upper-limit to 0 in problem
        # Store old limit for later to restore
        if partner_reaction is not None:
            partner_id = partner_reaction.id
            old_partner_ub = problem.variables.get_upper_bounds(partner_id)
            problem.variables.set_upper_bounds(partner_id, 0.0)

        utils.reset_objective(problem)
        problem.objective.set_linear(reaction.id, 1.0)
        problem.objective.set_name(str(reaction.id))
        problem.objective.set_sense(problem.objective.sense.maximize)

        problem.solve()
        rxn_max = problem.solution.get_objective_value()

        sub_cache[reaction.id] = rxn_max

        # Restore limit of partner reaction to old state
        if partner_reaction is not None:
            partner_id = partner_reaction.id
            problem.variables.set_upper_bounds(partner_id, old_partner_ub)

    return sub_cache

def maximize_metab_range(start_stop, args):
    """
    For precaching the maxmimum feasible exchange values of metabolites

    Returns
    -------
    sub_cache: dict
        key: species_id
        value: maximum flux
    """
    sub_cache = {}
    model = models.init_model(args['model'], species=args['species'],
                              exchange_limit=EXCHANGE_LIMIT,
                              media=args['media'])
    problem = initialize_cplex_problem(model, args['num_threads'], args['lpmethod'])

    metabolites = sorted(list(model.species.values()), key=lambda r:r.id)[start_stop[0]:start_stop[1]]


    all_names = set(problem.linear_constraints.get_names())

    for metabolite in tqdm(metabolites, file=sys.stderr):

        met_id = metabolite.id

        if met_id not in all_names:
            # This can happen if the metabolite does not participate
            # in any reaction. As a result, it won't be in any
            # constraints - happens in RECON2
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

        #if the selected_rxns or only_exchange options are used --> then we don't want to add reactions unless one of the pair already exists

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
        #secretion_max = maximize_reaction(model, problem, secretion_rxn)
        utils.reset_objective(problem)
        problem.objective.set_linear(secretion_rxn, 1.0)
        problem.objective.set_name(str(secretion_rxn))
        problem.objective.set_sense(problem.objective.sense.maximize)

        problem.solve()
        rxn_max = problem.solution.get_objective_value()

        sub_cache[secretion_rxn] = rxn_max

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
        #uptake_max = maximize_reaction(model, problem, uptake_rxn)
        utils.reset_objective(problem)
        problem.objective.set_linear(uptake_rxn, 1.0)
        problem.objective.set_name(str(uptake_rxn))
        problem.objective.set_sense(problem.objective.sense.maximize)

        problem.solve()
        rxn_max = problem.solution.get_objective_value()

        sub_cache[uptake_rxn] = rxn_max

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
            
    return sub_cache