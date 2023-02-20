from __future__ import print_function, division, absolute_import
import cplex
import argparse
import os
import multiprocessing
import numpy as np
import pandas as pd
import sys
import subprocess as sp
import logging
import datetime
import json
import gzip
import time
import timeit
from functools import partial
from tqdm import tqdm
from six import string_types
from math import ceil
from random import shuffle

from .algorithm import maximize_reaction, initialize_cplex_problem
from ..models import init_model

from .. import globals
from .. import utils
from .. import models
from . import cache

def runCompassParallelTransposed(args):
    logger = logging.getLogger('compass')

    # If we're here, then run compass on this machine with N processes
    if args['num_processes'] is None:
        args['num_processes'] = multiprocessing.cpu_count()

    if args['num_processes'] > multiprocessing.cpu_count():
        args['num_processes'] = multiprocessing.cpu_count()

    # Get the number of samples
    data = utils.read_data(args['data'])#pd.read_csv(args['data'], sep='\t', index_col=0)
    n_samples = len(data.columns)

    #   Get number of reactions and metabolites
    model = init_model(model=args['model'], species=args['species'],
                       exchange_limit=globals.EXCHANGE_LIMIT, media=args['media'], 
                       isoform_summing=args['isoform_summing'])
    n_reactions = len(model.reactions)
    if args['test_mode']:
        n_reactions = min(n_reactions, 100)
    n_metabolites = len(model.species)
    if args['test_mode']:
        n_metabolites = min(n_metabolites, 50)

    partial_map_fun = partial(_parallel_map_fun_transposed, args=args)

    pool = multiprocessing.Pool(args['num_processes'])
    sample_indices = np.linspace(0, n_samples, args['num_processes']+1, dtype=int) 
    sample_indices = [(sample_indices[i], sample_indices[i+1]) for i in range(len(sample_indices)-1)]
    ranges = [(i, [0, n_reactions], [0, n_metabolites], sample_indices[i]) for i in range(len(sample_indices))]

    logger.info(
        "Processing {} jobs using {} processes"
        .format(len(ranges), args['num_processes'])
    )

    logger.info(
        "Progress bar will update once the first job is finished"
    )

    logging_dir = os.path.join(args['temp_dir'], 'logging')
    if not os.path.isdir(logging_dir):
        os.makedirs(logging_dir)

    pbar = tqdm(total=len(ranges))

    for _ in pool.imap_unordered(partial_map_fun, ranges):
        pbar.update()

    collectCompassResultsTransposed(args['data'], args['temp_dir'],
                          args['output_dir'], args)

    logger.info("COMPASS Completed Successfully")

def _parallel_map_fun_transposed(ranges, args):
        i, ranges = ranges[0], ranges[1:]

        data = args['data']
        model = args['model']
        media = args['media']
        temp_dir = args['temp_dir']
        logging_dir = os.path.join(temp_dir, 'logging')

        out_file = os.path.join(logging_dir, 'out'+str(i)+'.log')
        err_file = os.path.join(logging_dir, 'err'+str(i)+'.log')
        with open(out_file, 'w') as fout, open(err_file, 'w') as ferr:
            stdout_bak = sys.stdout
            stderr_bak = sys.stderr
            sys.stdout = fout
            sys.stderr = ferr

            globals.init_logger(logging_dir)

            logger = logging.getLogger('compass')
            logger.debug("Compass: Single-sample mode")
            logger.debug("Supplied Arguments: ")
            for (key, val) in args.items():
                logger.debug("   {}: {}".format(key, val))
            logger.debug("Ranges: "+str(ranges))

            start_time = datetime.datetime.now()
            logger.debug("\nCOMPASS Started: {}".format(start_time))

            try:
                compass_transposed(ranges=ranges, 
                            data=data, model=model, 
                            media=media, args=args)
            except Exception as e:
                sys.stdout = stdout_bak
                sys.stderr = stderr_bak

                # Necessary because cplex exceptions can't be pickled
                #   and can't transfer from subprocess to main process
                if 'cplex' in str(type(e)).lower():
                    raise(Exception(str(e)))
                else:
                    raise(e)

            end_time = datetime.datetime.now()
            logger.debug("\nElapsed Time: {}".format(end_time-start_time))

def collectCompassResultsTransposed(data, temp_dir, out_dir, args):
    """
    Collects results for individual samples in temp_dir
    and aggregates into out_dir

    Parameters
    ==========
    data : str
       Full path to data file

    temp_dir : str
        Directory - where to look for sample results.

    out_dir : str
        Where to store aggregated results.  Is created if it doesn't exist.

    args : dict
        Other arguments
    """

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    logger = logging.getLogger('compass')
    logger.info("Collecting results from: " + temp_dir)
    logger.info("Writing output to: " + out_dir)

    # Get the number of samples
    expression = utils.read_data(data)#pd.read_csv(data, sep='\t', index_col=0)
    n_samples = len(expression.columns)

    reactions_all = []
    secretions_all = []
    uptake_all = []

    # Gather all the results
    for i in range(n_samples):

        sample_name = expression.columns[i]
        sample_dir = os.path.join(temp_dir, 'sample' + str(i))

        try:
            reactions = pd.concat([pd.read_csv(os.path.join(sample_dir, x), sep='\t', index_col=0) 
                for x in os.listdir(sample_dir) if 'reactions' == x[:9] and '.txt' == x[-4:]])

            reactions_all.append(reactions)

        except:
            reactions_all.append(pd.DataFrame(columns=[sample_name]))

        try:
            secretions = pd.concat([pd.read_csv(os.path.join(sample_dir, x), sep='\t', index_col=0) 
                for x in os.listdir(sample_dir) if 'secretions' == x[:9] and '.txt' == x[-4:]])

            secretions_all.append(secretions)

        except:
            secretions_all.append(pd.DataFrame(columns=[sample_name]))

        try:
            uptake = pd.concat([pd.read_csv(os.path.join(sample_dir, x), sep='\t', index_col=0) 
                for x in os.listdir(sample_dir) if 'uptake' == x[:9] and '.txt' == x[-4:]])

            uptake_all.append(uptake)

        except:
            uptake_all.append(pd.DataFrame(columns=[sample_name]))

    # Join and output
    if not args['no_reactions']:
        reactions_all = pd.concat(reactions_all, axis=1, sort=True)
        reactions_all.to_csv(
            os.path.join(out_dir, 'reactions.tsv'), sep="\t")

    if args['calc_metabolites']:
        secretions_all = pd.concat(secretions_all, axis=1, sort=True)
        secretions_all.to_csv(
            os.path.join(out_dir, 'secretions.tsv'), sep="\t")

        uptake_all = pd.concat(uptake_all, axis=1, sort=True)
        uptake_all.to_csv(
            os.path.join(out_dir, 'uptake.tsv'), sep="\t")

    # Output a JSON version of the model
    model = init_model(model=args['model'], species=args['species'],
                       exchange_limit=globals.EXCHANGE_LIMIT, media=args['media'], 
                       isoform_summing=args['isoform_summing'])

    model_file = os.path.join(out_dir, 'model.json.gz')

    with gzip.open(model_file, 'w') as gzfile:
        gzfile.write(model.to_JSON().encode('utf-8'))

def compass_transposed(ranges, data, model, media, args):
    logger = logging.getLogger("compass")

    model = init_model(model=args['model'], species=args['species'],
                       exchange_limit=globals.EXCHANGE_LIMIT, media=args['media'], 
                       isoform_summing=args['isoform_summing'])

    reaction_range, metabolite_range, sample_range = ranges
    if not args['reaction_range']:
        reaction_range = [0, len(model.getReactions())]
    if not args['metabolite_range']:
        metabolite_range = [0, model.species.values()]

    problem = initialize_cplex_problem(model, args['num_threads'], args['lpmethod'], args['advance'])

    expression = utils.read_data(data)
    if not args['sample_range']:
        sample_range = [0, len(expression.columns)]
    
    sample_names = [str(x) for x in expression.columns[sample_range[0]:sample_range[1]]]

    reaction_scores = {sample_name:{} for sample_name in sample_names}
    reaction_penalties = {}

    reaction_penalties = pd.read_csv(
            args['penalties_file'], sep="\t", header=0, index_col = 0)
    
    reactions = list(model.reactions.values())
    reaction_order = open('/data/yosef/users/schel337/reaction_order.txt', 'r').read().split()
    reaction_order = [int(x) for x in reaction_order]
    reactions = [reactions[x] for x in reaction_order]

    if args['test_mode']:
        reactions = reactions[0:100]
    reactions = reactions[reaction_range[0]:reaction_range[1]]

    directories = [os.path.join(args['temp_dir'], 'sample' + str(i)) for i in range(sample_range[0], sample_range[1])]

    for directory in directories:
        if not os.path.isdir(directory) and directory != '/dev/null':
            os.makedirs(directory)

        #TBD fix potential glitch, where different reactions from same sample are done
        if os.path.exists(os.path.join(directory, 'success_token')):
            logger.info('success_token detected, results already calculated.')
            logger.info('COMPASS Completed Successfully')
            return

    model_cache = cache.load(model)

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
        
        r_max = maximize_reaction(model, problem, reaction.id, use_cache = True)

        if r_max == 0:
            for sample_name in sample_names:
                reaction_scores[sample_name][reaction.id] = 0
            # If Reaction can't carry flux anyways, just continue
            
        else:
            problem.linear_constraints.add(
                    lin_expr=[cplex.SparsePair(ind=[reaction.id], val=[1.0])],
                    senses=['R'],
                    rhs=[globals.BETA * r_max],
                    names=['REACTION_OPT'])

            for sample_name in sample_names:
                # Minimize Penalty
                utils.reset_objective(problem)
                problem.objective.set_linear(
                    list(reaction_penalties[sample_name].iteritems())
                )
                problem.objective.set_sense(problem.objective.sense.minimize)
                
                problem.solve()

                value = problem.solution.get_objective_value()
                reaction_scores[sample_name][reaction.id] = value

            # Remove Constraint
            problem.linear_constraints.delete('REACTION_OPT')

        if partner_reaction is not None:
            partner_id = partner_reaction.id
            problem.variables.set_upper_bounds(partner_id, old_partner_ub)

    # Output results to file
    logger.info("Writing output files...")
    if not args['no_reactions']:
        for i in range(len(sample_names)):
            reaction_scores[sample_names[i]] = pd.Series(reaction_scores[sample_names[i]], name=sample_names[i]).sort_index()
            reaction_scores[sample_names[i]].to_csv(os.path.join(directories[i], 'reactions'+str(reaction_range)+'.txt'),
                                sep="\t", header=True)

    # write success token. Doesn't work when sampels are done with multiple calls
    #with open(os.path.join(directory, 'success_token'), 'w') as fout:
    #    fout.write('Success!')
    
    logger.info('COMPASS Completed Successfully')