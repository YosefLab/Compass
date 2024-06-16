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
from ..globals import BETA, EXCHANGE_LIMIT, LICENSE_DIR
import compass.global_state as global_state

import gurobipy as gp
from gurobipy import GRB

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

    model = models.init_model(model=args['model'], species=args['species'],
                       exchange_limit=EXCHANGE_LIMIT, media=args['media'], 
                       isoform_summing=args['isoform_summing'])

    logger.info("Running COMPASS on model: %s", model.name)

    perf_log = None
    if args['detailed_perf']:
        cols = ['order','max rxn time', 'max rxn method', 'cached', 'min penalty time', 
        'min penalty method', 'min penalty sensitvivity', 'kappa']
        perf_log = {c:{} for c in cols}

    if args['generate_cache']:
        cache.clear(model) #TBD add media specifier here too

    # Build model into Gurobi model
    credentials = utils.parse_gurobi_license_file(os.path.join(LICENSE_DIR, 'gurobi.lic'))
    gp_model = initialize_gurobi_model(model, credentials, args['num_threads'], args['lpmethod'], args['advance'])

    # Only read this to get the number of samples and the sample name
    # Use nrows=1 so this is fast
    samples = utils.read_sample_names(data)
    if samples is None:
        sample_name = 'sample_'+str(sample_index)
        logger.info("Processing Sample %i: %s", sample_index, sample_name)
    else:
        sample_name = str(samples[sample_index])
        logger.info("Processing Sample %i/%i: %s", sample_index,
            len(samples), sample_name)
    global_state.set_current_cell_name(sample_name)

    # Run core compass algorithm

    # Evaluate reaction penalties
    # Actual time spent on evaluating reaction penalties is not calculated
    # Logger actually records time used to read already computed reaction penalties
    penalty_start = time.process_time()
    logger.info("Evaluating Reaction Penalties...")
    reaction_penalties = pd.read_csv(
        args['penalties_file'], sep="\t", header=0,
        usecols=[0, sample_index + 1]) #0 is the Reaction column,

    reaction_penalties = reaction_penalties.set_index("Reaction").iloc[:, 0]
    penalty_elapsed = time.process_time() - penalty_start

    react_start = time.process_time()
    if not args['no_reactions']:
        logger.info("Evaluating Reaction Scores...")
        reaction_scores = compass_reactions(
            model, gp_model, reaction_penalties,
            perf_log=perf_log, args=args)
    react_elapsed = time.process_time() - react_start

    #if user wants to calc reaction scores, but doesn't want to calc metabolite scores, calc only the exchange reactions
    logger.info("Evaluating Exchange/Secretion/Uptake Scores...")
    exchange_start = time.process_time()
    uptake_scores, secretion_scores, exchange_rxns = compass_exchange(
        model, gp_model, reaction_penalties,
        only_exchange=(not args['no_reactions']) and not args['calc_metabolites'],
        perf_log=perf_log, args=args)
    exchange_elapsed = time.process_time() - exchange_start

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

    logger.info("Compass Penalty Time: "+str(penalty_elapsed))
    if not args['no_reactions']:
        logger.info("Compass Reaction Time: "+str(react_elapsed))
        logger.info("Processed "+str(len(reaction_scores))+" reactions")
    logger.info("Compass Exchange Time: "+str(exchange_elapsed))
    logger.info("Processed "+str(len(uptake_scores))+" uptake reactions")
    logger.info("Processed "+str(len(secretion_scores))+" secretion reactions")
    
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


def compass_exchange(model, gp_model, reaction_penalties, only_exchange=False, perf_log=None, args = None):
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

    # Setting only_exchange=False might cause a pair of uptake/secretion reactions associated with a certain
    # metabolite to be added even if there was no uptake/secretion reaction associated with this metabolite
    # only_exchange=True adds the other half of the pair if there is an uptake/secretion reaction and does not
    # if there is no such reaction

    secretion_scores = {}
    uptake_scores = {}
    exchange_rxns = {}
    metabolites = list(model.species.values())
    if args['test_mode']:
        metabolites = metabolites[0:50]

    #populate the list of selected_reaction_ids - do this once outside of the loop
    if args['select_reactions'] or args['select_subsystems']:
        selected_reaction_ids = read_selected_reactions(args['select_reactions'], args['select_subsystems'], model)

    for metabolite in tqdm(metabolites, file=sys.stderr):

        met_id = metabolite.id
        met_id_constr = gp_model.getConstrByName(met_id)

        if met_id_constr is None:
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
        sp = gp_model.getRow(met_id_constr)
        rxn_ids = [sp.getVar(i).varname for i in range(sp.size())]
        reactions = [model.reactions[x] for x in rxn_ids]

        #If user wants only exchange reaction - limit the reactions space through which we iterate
        if only_exchange:
            reactions = [x for x in reactions if x.is_exchange]

        if args['select_reactions'] or args['select_subsystems']:
            #r.id is a unidirectional identifier (ending with _pos or _neg suffix --> we remove it and compare to the undirected reaction id)
            reactions = [r for r in reactions if ((r.id)[:-4] in selected_reaction_ids or str(r.id) in selected_reaction_ids)]

        # Extra reactions are duplicates
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
            rxn_index = gp_model.addVar(
                name=secretion_rxn,
                ub=model.maximum_flux,
                lb=0.0,
                vtype=GRB.CONTINUOUS,
            )

            # Add it to the metabolites constraint
            gp_model.chgCoeff(met_id_constr, rxn_index, -1.0)
            gp_model.update()

        #if only exchange flag is set - don't add uptakes that do not exist
        if (uptake_rxn is None):
            added_uptake = True
            uptake_rxn = met_id + "_UPTAKE"

            # Add uptake reaction to the problem as a variable
            rxn_index = gp_model.addVar(
                name=uptake_rxn,
                ub=EXCHANGE_LIMIT,
                lb=0.0,
                vtype=GRB.CONTINUOUS,
            )

            # Add it to the metabolite's constraint
            # TODO: need to check if met_id_constr is modified by secrete_rxn
            # Find way to remove gp_model.update() in secrete_rxn
            gp_model.chgCoeff(met_id_constr, rxn_index, 1.0)
            gp_model.update()

        # Modify the constraint in the problem
        #   e.g. Add the metabolites connections
        all_uptake = [uptake_rxn] + extra_uptake_rxns
        all_secretion = [secretion_rxn] + extra_secretion_rxns

        # -----------------
        # Optimal Secretion
        # -----------------

        # Close all uptake, storing their upper-bounds to restore later
        old_uptake_upper = {}
        for rxn_id in all_uptake:
            rxn_var = gp_model.getVarByName(rxn_id)
            old_ub = rxn_var.ub
            old_uptake_upper[rxn_id] = old_ub
            old_lb = rxn_var.lb
            rxn_var.setAttr('ub', max(old_lb, 0))

        # Close extra secretion, storing upper-bounds to restore later
        old_secretion_upper = {}
        for rxn_id in extra_secretion_rxns:
            rxn_var = gp_model.getVarByName(rxn_id)
            old_ub = rxn_var.ub
            old_secretion_upper[rxn_id] = old_ub
            old_lb = rxn_var.lb
            rxn_var.setAttr('ub', max(old_lb, 0))

        # Get max of secretion reaction
        secretion_max = maximize_reaction(model, gp_model, secretion_rxn, perf_log=perf_log)

        # Set contraint of max secretion to BETA*max
        secretion_var = gp_model.getVarByName(secretion_rxn)
        rhs = BETA * secretion_max
        name = 'SECRETION_OPT'
        gp_model.addConstr(secretion_var >= rhs, name=name)

        # Find minimimum penalty
        obj = gp.LinExpr()
        # list contains elements of form (reaction, penalty) associated with the current sample
        # This sets the objective to be v_1 * p_1 + v_2 * p_2 + ... + v_m * p_m
        for rxn, penalty in reaction_penalties.items():
            obj += penalty * gp_model.getVarByName(rxn)
        gp_model.setObjective(obj, GRB.MINIMIZE)
        gp_model.update()
        
        if perf_log is not None:
            start_time = time.process_time()
        #gp_model.optimize()
        if perf_log is not None:
            perf_log['min penalty time'][secretion_rxn] = time.process_time() - start_time
        #value = gp_model.ObjVal
        global_state.set_current_reaction_id(secretion_rxn)
        value = optimize_model_wrapper(gp_model)
        secretion_scores[met_id] = value

        # Clear Secretion constraint
        constr = gp_model.getConstrByName('SECRETION_OPT')
        gp_model.remove(constr)

        # Restore all uptake
        for rxn_id, old_ub in old_uptake_upper.items():
            uptake_var = gp_model.getVarByName(rxn_id)
            uptake_var.setAttr('ub', old_ub)

        # Restore extra secretion
        for rxn_id, old_ub in old_secretion_upper.items():
            secretion_var = gp_model.getVarByName(rxn_id)
            secretion_var.setAttr('ub', old_ub)

        # -----------------
        # Optimal Uptake
        # -----------------

        # Close extra uptake
        old_uptake_upper = {}
        for rxn_id in extra_uptake_rxns:
            rxn_var = gp_model.getVarByName(rxn_id)
            old_ub = rxn_var.ub
            old_uptake_upper[rxn_id] = old_ub
            old_lb = rxn_var.lb
            rxn_var.setAttr('ub', max(old_lb, 0))

        # Close all secretion
        old_secretion_upper = {}
        for rxn_id in all_secretion:
            rxn_var = gp_model.getVarByName(rxn_id)
            old_ub = rxn_var.ub
            old_secretion_upper[rxn_id] = old_ub
            old_lb = rxn_var.lb
            rxn_var.setAttr('ub', max(old_lb, 0))

        # Get max of uptake reaction
        uptake_max = maximize_reaction(model, gp_model, uptake_rxn, perf_log=perf_log)

        # Set contraint of max uptake with BETA*max
        uptake_var = gp_model.getVarByName(uptake_rxn)
        rhs = BETA * uptake_max
        name = 'UPTAKE_OPT'
        gp_model.addConstr(uptake_var >= rhs, name=name)

        # Find minimimum penalty
        obj = gp.LinExpr()
        # list contains elements of form (reaction, penalty) associated with the current sample
        # This sets the objective to be v_1 * p_1 + v_2 * p_2 + ... + v_m * p_m
        for rxn, penalty in reaction_penalties.items():
            obj += penalty * gp_model.getVarByName(rxn)
        gp_model.setObjective(obj, GRB.MINIMIZE)
        gp_model.update()

        if perf_log is not None:
            start_time = time.process_time()
        #gp_model.optimize()
        if perf_log is not None:
            perf_log['min penalty time'][uptake_rxn] = time.process_time() - start_time
        #value = gp_model.ObjVal
        global_state.set_current_reaction_id(uptake_rxn)
        value = optimize_model_wrapper(gp_model)
        uptake_scores[met_id] = value

        # Clear Secretion constraint
        constr = gp_model.getConstrByName('UPTAKE_OPT')
        gp_model.remove(constr)

        # Restore extra uptake
        for rxn_id, old_ub in old_uptake_upper.items():
            uptake_var = gp_model.getVarByName(rxn_id)
            uptake_var.setAttr('ub', old_ub)

        # Restore all secretion
        for rxn_id, old_ub in old_secretion_upper.items():
            secretion_var = gp_model.getVarByName(rxn_id)
            secretion_var.setAttr('ub', old_ub)

        # Remove added uptake and secretion reactions
        if added_uptake:
            uptake_var = gp_model.getVarByName(uptake_rxn)
            gp_model.remove(uptake_var)
        else:
            for rxn_id in all_uptake:
                exchange_rxns[rxn_id] = uptake_scores[met_id]

        if added_secretion:
            secretion_var = gp_model.getVarByName(secretion_rxn)
            gp_model.remove(secretion_var)
        else:
            for rxn_id in all_secretion:
                exchange_rxns[rxn_id] = secretion_scores[met_id]

    return uptake_scores, secretion_scores, exchange_rxns



def compass_reactions(model, gp_model, reaction_penalties, perf_log=None, args = None):

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
            partner_var = gp_model.getVarByName(partner_id)
            old_partner_ub = partner_var.ub
            old_partner_lb = partner_var.lb
            partner_var.setAttr('ub', max(old_partner_lb, 0))
        
        r_max = maximize_reaction(model, gp_model, reaction.id, perf_log=perf_log)

        # If Reaction can't carry flux anyways (v_r^opt = 0), just continue
        if r_max == 0:
            reaction_scores[reaction.id] = 0
            if perf_log is not None:
               perf_log['min penalty time'][reaction.id] = 0
               #perf_log['blocked'][reaction.id] = True

        else:
            expr = gp.LinExpr()
            expr += gp_model.getVarByName(reaction.id)
            rhs = BETA * r_max
            name = 'REACTION_OPT'
            gp_model.addConstr(expr >= rhs, name=name)

            # Minimize Penalty
            # TODO: does not check if objective name is 'reaction_penalties'
            # because gurobi does not assign names to objectives
            obj = gp.LinExpr()
            # list contains elements of form (reaction, penalty) associated with the current sample
            # This sets the objective to be v_1 * p_1 + v_2 * p_2 + ... + v_m * p_m
            for rxn, penalty in reaction_penalties.items():
                obj += penalty * gp_model.getVarByName(rxn)
            gp_model.setObjective(obj, GRB.MINIMIZE)
            gp_model.update()

            if perf_log is not None:
                #perf_log['blocked'][reaction.id] = False
                start_time = time.process_time()

            #gp_model.optimize()

            # TODO: modify for gurobi
            if perf_log is not None:
                perf_log['min penalty time'][reaction.id] = time.process_time() - start_time
                perf_log['min penalty method'][reaction.id] = gp_model.getParamInfo('Method')
                #perf_log['min penalty sensitvivity'][reaction.id] = problem.solution.sensitivity.objective(reaction.id)
                #if hasattr(problem.solution.get_quality_metrics(),'kappa'):
                   #perf_log['kappa'][reaction.id] = problem.solution.get_quality_metrics().kappa

            if args['save_argmaxes']:
                all_vars = gp_model.getVars()
                values = gp_model.getAttr("X", all_vars)
                argmaxes.append(np.array(values))
                argmaxes_order.append(reaction.id)

            #value = gp_model.ObjVal
            global_state.set_current_reaction_id(reaction.id)
            value = optimize_model_wrapper(gp_model)
            reaction_scores[reaction.id] = value

            # Remove Constraint
            constr = gp_model.getConstrByName('REACTION_OPT')
            gp_model.remove(constr)

        # Restore limit of partner reaction to old state
        if partner_reaction is not None:
            partner_id = partner_reaction.id

            partner_var = gp_model.getVarByName(partner_id)
            partner_var.setAttr('ub', old_partner_ub)

    if args['save_argmaxes']:
        argmaxes = np.vstack(argmaxes)
        np.save(os.path.join(args['save_argmaxes_dir'],'argmaxes.npy'), argmaxes)
        argmaxes_order = np.array(argmaxes_order)
        np.save(os.path.join(args['save_argmaxes_dir'],'argmaxes_order.npy'), argmaxes_order)

    return reaction_scores


def initialize_gurobi_model(model, credentials, num_threads=1, lpmethod=-1, adv=2):
    # type: (compass.models.MetabolicModel)
    """
    Builds and returns a gurobi model representing our metabolic model

    Limits exchange reactions and makes all reactions unidirectional
    by splitting into two components
    """

    # Create the Gurobi model
    env = gp.Env(params=credentials)
    gp_model = gp.Model(env=env)

    # Set Parameters for the Gurobi model

    gp_model.setParam("OutputFlag", 0)  # Disable all output
    gp_model.setParam("LogToConsole", 0)  # Disable console output

    # Set numerical emphasis to improve precision
    gp_model.setParam("NumericFocus", 3)  # Equivalent to numerical emphasis in CPLEX

    # Set number of threads
    gp_model.setParam("Threads", num_threads)  # Set the number of threads to use

    # Set the primal and dual preprocessing options
    gp_model.setParam("Presolve", 2)  # 2 means aggressive presolve, 1 for conservative, 0 for off

    # Set optimization method
    gp_model.setParam("Method", lpmethod)  # 0: Automatic, 1: Primal Simplex, 2: Dual Simplex, etc.

    # Set optimality tolerance
    gp_model.setParam("OptimalityTol", 1e-9)  # Default is 1e-6, minimum is 1e-9

    # Set barrier convergence tolerance
    gp_model.setParam("BarConvTol", 1e-12)  # Default is 1e-8, minimum is 1e-12

    gp_model.setParam(GRB.Param.Threads, num_threads)
    gp_model.setParam(GRB.Param.Method, lpmethod)

    # Add variables
    reactions = list(model.reactions.values())

    # Define minimum and maximum flux for each reaction
    for x in reactions:
        gp_model.addVar(
            lb=x.lower_bound, 
            ub=x.upper_bound, 
            name=x.id, 
            vtype=GRB.CONTINUOUS)
    gp_model.update()

    # Add constraints

    # Add stoichiometry constraints

    '''
    utils.get_steadystate_constraints

    For each metabolite:
        c_lin_expr is the zero flux linear expression, of form
            c_1 * r_1 + c_2 * r_2 + ... + c_m * r_m = 0
        sense is 'E', representing equal to 0
        rhs is 0
        name of corresponding metabolite is given to entire linear expression, not variables
    '''

    c_lin_expr, c_rhs, c_names = (
        utils.get_steadystate_constraints(model, gp_model))

    for lin_expr, rhs, name in zip(c_lin_expr, c_rhs, c_names):
        gp_model.addConstr(lin_expr == rhs, name=name)
    gp_model.update()

    # Initialize the objective
    ### utils.reset_objective(gp_model)

    return gp_model


def maximize_reaction(model, gp_model, rxn, use_cache=True, perf_log=None):

    """Maximizes the current reaction in the problem
    Attempts to retrieve the value from cache if its in cache
    """
    
    if perf_log is not None:
        start_time = time.process_time()
        perf_log['order'][rxn] = len(perf_log['order'])

    # Load from cache if it exists and return
    if use_cache:
        model_cache = cache.load(model)
        if rxn in model_cache:
            if perf_log is not None:
                perf_log['cached'][rxn] = True
                perf_log['max rxn time'][rxn] = time.process_time() - start_time
            return model_cache[rxn]

    # Maximize the reaction
    ### utils.reset_objective(gp_model)
    gp_model.setObjective(gp_model.getVarByName(rxn), GRB.MAXIMIZE)
    gp_model.update()

    #gp_model.optimize()
    #rxn_max = gp_model.ObjVal
    global_state.set_current_reaction_id(rxn)
    rxn_max = optimize_model_wrapper(gp_model)

    # Save the result
    model_cache = cache.load(model)
    model_cache[rxn] = rxn_max

    if perf_log is not None:
        perf_log['cached'][rxn] = False
        perf_log['max rxn time'][rxn] = time.process_time() - start_time
        # TODO: modify for gurobi
        #perf_log['max rxn method'][rxn] = problem.solution.get_method()
        perf_log['max rxn method'][rxn] = gp_model.getParamInfo('Method')

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
    model = models.init_model(model=args['model'], species=args['species'],
                       exchange_limit=EXCHANGE_LIMIT, media=args['media'], 
                       isoform_summing=args['isoform_summing'])
    credentials = utils.parse_gurobi_license_file(os.path.join(LICENSE_DIR, 'gurobi.lic'))
    gp_model = initialize_gurobi_model(model, credentials, args['num_threads'], args['lpmethod'])

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
            partner_var = gp_model.getVarByName(partner_id)
            old_partner_ub = partner_var.ub
            old_partner_lb= partner_var.lb
            partner_var.setAttr('ub', max(old_partner_lb, 0))
        

        utils.reset_objective(gp_model)
        gp_model.setObjective(gp_model.getVarByName(reaction.id), GRB.MAXIMIZE)
        gp_model.update()

        gp_model.optimize()
        rxn_max = gp_model.ObjVal

        sub_cache[reaction.id] = rxn_max

        # Restore limit of partner reaction to old state
        if partner_reaction is not None:
            partner_id = partner_reaction.id

            partner_var = gp_model.getVarByName(partner_id)
            partner_var.setAttr('ub', old_partner_ub)

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
    model = models.init_model(model=args['model'], species=args['species'],
                       exchange_limit=EXCHANGE_LIMIT, media=args['media'], 
                       isoform_summing=args['isoform_summing'])
    credentials = utils.parse_gurobi_license_file(os.path.join(LICENSE_DIR, 'gurobi.lic'))
    gp_model = initialize_gurobi_model(model, credentials, args['num_threads'], args['lpmethod'])

    # init_model returns entire list of metabolites as specified by the RECON2 model
    # However, some metabolites are not associated with any reactions,
    # and maximize_metab_range only computes those that are associated with reactions
    # Therefore model.species might contain metabolites that don't need to be maximized at all
    metabolites = sorted(list(model.species.values()), key=lambda r:r.id)[start_stop[0]:start_stop[1]]

    for metabolite in tqdm(metabolites, file=sys.stderr):

        met_id = metabolite.id
        met_id_constr = gp_model.getConstrByName(met_id)

        if met_id_constr is None:
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
        sp = gp_model.getRow(met_id_constr)
        rxn_ids = [sp.getVar(i).varname for i in range(sp.size())]
        reactions = [model.reactions[x] for x in rxn_ids]

        # NOTE: only_exchange flag does not apply to maximize_metab_range
        # however this has no effect since filtering for exchange reactions
        # is still done within 'for reaction in reactions' loop

        # if only_exchange:
        #     reactions = [x for x in reactions if x.is_exchange]

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

        # NOTE: only_exchange flag does not apply to maximize_metab_range
        # i.e. pairs of exchange reactions can be added even if neither exist

        #if the selected_rxns or only_exchange options are used --> then we don't want to add reactions unless one of the pair already exists

        if (secretion_rxn is None):
            added_secretion = True
            secretion_rxn = met_id + "_SECRETION"

            # Add secretion reaction to the problem as a variable
            rxn_index = gp_model.addVar(
                name=secretion_rxn,
                ub=model.maximum_flux,
                lb=0.0,
                vtype=GRB.CONTINUOUS,
            )

            # Add it to the metabolites constraint
            gp_model.chgCoeff(met_id_constr, rxn_index, -1.0)
            gp_model.update()


        #if only exchange flag is set - don't add uptakes that do not exist
        if (uptake_rxn is None):
            added_uptake = True
            uptake_rxn = met_id + "_UPTAKE"

            # Add uptake reaction to the problem as a variable
            rxn_index = gp_model.addVar(
                name=uptake_rxn,
                ub=EXCHANGE_LIMIT,
                lb=0.0,
                vtype=GRB.CONTINUOUS,
            )

            # Add it to the metabolite's constraint
            # TODO: need to check if met_id_constr is modified by secrete_rxn
            # Find way to remove gp_model.update() in secrete_rxn
            gp_model.chgCoeff(met_id_constr, rxn_index, 1.0)
            gp_model.update()

        # Modify the constraint in the problem
        #   e.g. Add the metabolites connections
        all_uptake = [uptake_rxn] + extra_uptake_rxns
        all_secretion = [secretion_rxn] + extra_secretion_rxns

        # -----------------
        # Optimal Secretion
        # -----------------

        # Close all uptake, storing their upper-bounds to restore later
        old_uptake_upper = {}
        for rxn_id in all_uptake:
            rxn_var = gp_model.getVarByName(rxn_id)
            old_ub = rxn_var.ub
            old_uptake_upper[rxn_id] = old_ub
            old_lb = rxn_var.lb
            rxn_var.setAttr('ub', max(old_lb, 0))

        # Close extra secretion, storing upper-bounds to restore later
        old_secretion_upper = {}
        for rxn_id in extra_secretion_rxns:
            rxn_var = gp_model.getVarByName(rxn_id)
            old_ub = rxn_var.ub
            old_secretion_upper[rxn_id] = old_ub
            old_lb = rxn_var.lb
            rxn_var.setAttr('ub', max(old_lb, 0))

        # Get max of secretion reaction
        #secretion_max = maximize_reaction(model, problem, secretion_rxn)
        utils.reset_objective(gp_model)
        gp_model.setObjective(gp_model.getVarByName(secretion_rxn), GRB.MAXIMIZE)

        gp_model.optimize()
        rxn_max = gp_model.ObjVal

        # NOTE: since only_exchange flag is not applicable here, additional secretion reactions
        # can be created and added to the cache
        sub_cache[secretion_rxn] = rxn_max

        # Restore all uptake
        for rxn_id, old_ub in old_uptake_upper.items():
            uptake_var = gp_model.getVarByName(rxn_id)
            uptake_var.setAttr('ub', old_ub)

        # Restore extra secretion
        for rxn_id, old_ub in old_secretion_upper.items():
            secretion_var = gp_model.getVarByName(rxn_id)
            secretion_var.setAttr('ub', old_ub)

        # -----------------
        # Optimal Uptake
        # -----------------

        # Close extra uptake
        old_uptake_upper = {}
        for rxn_id in extra_uptake_rxns:
            rxn_var = gp_model.getVarByName(rxn_id)
            old_ub = rxn_var.ub
            old_uptake_upper[rxn_id] = old_ub
            old_lb = rxn_var.lb
            rxn_var.setAttr('ub', max(old_lb, 0))

        # Close all secretion
        old_secretion_upper = {}
        for rxn_id in all_secretion:
            rxn_var = gp_model.getVarByName(rxn_id)
            old_ub = rxn_var.ub
            old_secretion_upper[rxn_id] = old_ub
            old_lb = rxn_var.lb
            rxn_var.setAttr('ub', max(old_lb, 0))

        # Get max of uptake reaction
        #uptake_max = maximize_reaction(model, problem, uptake_rxn)
        utils.reset_objective(gp_model)
        gp_model.setObjective(gp_model.getVarByName(uptake_rxn), GRB.MAXIMIZE)

        gp_model.optimize()
        rxn_max = gp_model.ObjVal

        # NOTE: since only_exchange flag is not applicable here, additional uptake reactions
        # can be created and added to the cache
        sub_cache[uptake_rxn] = rxn_max

        # Restore extra uptake
        for rxn_id, old_ub in old_uptake_upper.items():
            uptake_var = gp_model.getVarByName(rxn_id)
            uptake_var.setAttr('ub', old_ub)

        # Restore all secretion
        for rxn_id, old_ub in old_secretion_upper.items():
            secretion_var = gp_model.getVarByName(rxn_id)
            secretion_var.setAttr('ub', old_ub)

        # Remove added uptake and secretion reactions
        if added_uptake:
            uptake_var = gp_model.getVarByName(uptake_rxn)
            gp_model.remove(uptake_var)

        if added_secretion:
            secretion_var = gp_model.getVarByName(secretion_rxn)
            gp_model.remove(secretion_var)
            
    return sub_cache


def optimize_model_wrapper(gp_model) -> float:
    r"""
    Only optimize the Gurobi model if the reaction is selected for the cell. Else,
    skip the computation and return np.nan.
    """
    if global_state.current_reaction_is_selected_for_current_cell():
        gp_model.optimize()
        return gp_model.ObjVal
    else:
        return np.nan