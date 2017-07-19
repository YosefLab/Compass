"""
Entry points for compass
"""
from __future__ import absolute_import, print_function, division
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
from functools import partial
from tqdm import tqdm

from .._version import __version__
from .torque import submitCompassTorque
from .algorithm import singleSampleCompass
from .. import globals


def parseArgs():
    """Defines the command-line arguments and parses the Compass call

    Returns
    -------
    argparse.Namespace

    """
    parser = argparse.ArgumentParser(
                        prog="Compass",
                        description="Metabolic Modeling for Single Cells")

    parser.add_argument("--data", help="Gene expression matrix",
                        required=True,
                        metavar="FILE")

    parser.add_argument("--model", help="Metabolic Model to Use",
                        default="RECON2_mat",
                        choices=["RECON1_mat", "RECON2_mat"],
                        metavar="MODEL")

    parser.add_argument("--species",
                        help="Species to use to match genes to model",
                        choices=["homo_sapiens", "mus_musculus"],
                        metavar="SPECIES",
                        default="homo_sapiens")

    parser.add_argument("--media", help="Which media to simulate",
                        metavar="MEDIA")

    parser.add_argument("--output-dir", help="Where to store outputs",
                        default='.',
                        metavar="DIR")

    parser.add_argument("--temp-dir", help="Where to store temporary files",
                        default='<output-dir>/_tmp',
                        metavar="DIR")

    parser.add_argument("--torque-queue", help="Submit to a Torque queue",
                        metavar="QUEUE")

    parser.add_argument("--num-processes",
                        help="Limit to <N> Processes.  "
                             "Ignored when submitting job onto a queue",
                        type=int,
                        metavar="N")

    parser.add_argument("--lambda",
                        help="Smoothing factor for single-cell data. Should be"
                        " set between 0 and 1",
                        type=float,
                        default=0,
                        metavar="F")

    parser.add_argument("--single-sample",
                        help=argparse.SUPPRESS,
                        type=int,
                        metavar="N")

    parser.add_argument(
        "--and-function",
        help="Which function used to aggregate AND associations",
        choices=["min", "median", "mean"],
        metavar="FXN",
        default="min")

    # Hidden argument.  Used for batch jobs
    parser.add_argument("--collect", action="store_true",
                        help=argparse.SUPPRESS)

    parser.add_argument("--perplexity",
                        help="Effective number of neighbors for tsne kernel",
                        default=30,
                        type=int,
                        metavar="N")

    parser.add_argument("--symmetric-kernel", action="store_true",
                        help="Use symmetric TSNE kernel (slower)")

    # Also used for batch jobs
    parser.add_argument("--config-file", help=argparse.SUPPRESS)

    args = parser.parse_args()

    args = vars(args)  # Convert to a Dictionary

    load_config(args)

    # Convert directories/files to absolute paths
    args['data'] = os.path.abspath(args['data'])

    if args['temp_dir'] == "<output-dir>/_tmp":
        args['temp_dir'] = os.path.join(args['output_dir'], '_tmp')

    args['output_dir'] = os.path.abspath(args['output_dir'])
    args['temp_dir'] = os.path.abspath(args['temp_dir'])

    if args['lambda'] < 0 or args['lambda'] > 1:
        parser.error(
            "'lambda' parameter cannot be less than 0 or greater than 1"
        )

    return args


def entry():
    """Entry point for the compass command-line script
    """
    start_time = datetime.datetime.now()

    args = parseArgs()

    if not os.path.isdir(args['output_dir']):
        os.makedirs(args['output_dir'])

    # Log some things for debugging/record
    globals.init_logger(args['output_dir'])
    logger = logging.getLogger('mflux')
    logger.debug("MFlux version: " + __version__)

    try:
        commit = sp.check_output(
            ["git", '--git-dir', globals.GIT_DIR, "rev-parse", "--short",
             "HEAD"],
            stderr=open(os.devnull, 'w')
        )
        logger.debug("Git commit: " + commit.decode())
    except sp.CalledProcessError:
        logger.debug("Git commit: Not in Git repo")

    logger.debug("Python Version:")
    logger.debug(sys.version)
    logger.debug("Python prefix: " + sys.prefix)
    logger.debug("Numpy version: " + np.__version__)
    logger.debug("Pandas version: " + pd.__version__)
    logger.debug("Supplied Arguments: ")
    for (key, val) in args.items():
        logger.debug("   {}: {}".format(key, val))

    logger.debug("\nCOMPASS Started: {}".format(start_time))
    # Parse arguments and decide what course of action to take

    if args['single_sample'] is not None:
        singleSampleCompass(data=args['data'], model=args['model'],
                            media=args['media'], directory=args['temp_dir'],
                            sample_index=args['single_sample'], args=args)
        end_time = datetime.datetime.now()
        logger.debug("\nElapsed Time: {}".format(end_time-start_time))
        return

    if args['torque_queue'] is not None:
        submitCompassTorque(args,
                            temp_dir=args['temp_dir'],
                            output_dir=args['output_dir'],
                            queue=args['torque_queue'])
        return

    if args['collect']:
        collectCompassResults(args['data'], args['temp_dir'],
                              args['output_dir'])
        end_time = datetime.datetime.now()
        logger.debug("\nElapsed Time: {}".format(end_time-start_time))
        return

    # If we're here, then run compass on this machine with N processes
    if args['num_processes'] is None:
        args['num_processes'] = multiprocessing.cpu_count()

    if args['num_processes'] > multiprocessing.cpu_count():
        args['num_processes'] = multiprocessing.cpu_count()

    # Get the number of samples
    data = pd.read_table(args['data'], index_col=0)
    n_samples = len(data.columns)

    partial_map_fun = partial(_parallel_map_fun, args=args)

    pool = multiprocessing.Pool(args['num_processes'])

    logger.info(
        "Processing {} samples using {} processes"
        .format(n_samples, args['num_processes'])
    )

    logger.info(
        "Progress bar will update once the first sample is finished"
    )

    pbar = tqdm(total=n_samples)

    for _ in pool.imap_unordered(partial_map_fun, range(n_samples)):
        pbar.update()

    collectCompassResults(args['data'], args['temp_dir'], args['output_dir'])

    end_time = datetime.datetime.now()
    logger.debug("\nCompleted At: {}".format(end_time))
    logger.debug("\nElapsed Time: {}".format(end_time-start_time))
    logger.info(
        "COMPASS Completed Successfully"
    )


def _parallel_map_fun(i, args):

        data = args['data']
        model = args['model']
        media = args['media']
        temp_dir = args['temp_dir']

        sample_dir = os.path.join(temp_dir, 'sample' + str(i))

        if not os.path.isdir(sample_dir):
            os.makedirs(sample_dir)

        out_file = os.path.join(sample_dir, 'out.log')
        err_file = os.path.join(sample_dir, 'err.log')
        with open(out_file, 'w') as fout, open(err_file, 'w') as ferr:
            sys.stdout = fout
            sys.stderr = ferr

            globals.init_logger(sample_dir)

            logger = logging.getLogger('mflux')
            logger.debug("MFlux: Single-sample mode")
            logger.debug("Supplied Arguments: ")
            for (key, val) in args.items():
                logger.debug("   {}: {}".format(key, val))

            start_time = datetime.datetime.now()
            logger.debug("\nCOMPASS Started: {}".format(start_time))

            singleSampleCompass(
                data=data, model=model,
                media=media, directory=sample_dir,
                sample_index=i, args=args
            )

            end_time = datetime.datetime.now()
            logger.debug("\nElapsed Time: {}".format(end_time-start_time))


def collectCompassResults(data, temp_dir, out_dir):
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
    """

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    logger = logging.getLogger('mflux')
    logger.info("Collecting results from: " + temp_dir)
    logger.info("Writing output to: " + out_dir)

    # Get the number of samples
    expression = pd.read_table(data, index_col=0)
    n_samples = len(expression.columns)

    reactions_all = []
    secretions_all = []
    uptake_all = []

    # Gather all the results
    for i in range(n_samples):

        sample_name = expression.columns[i]
        sample_dir = os.path.join(temp_dir, 'sample' + str(i))

        try:
            reactions = pd.read_table(
                os.path.join(sample_dir, 'reactions.txt'),
                index_col=0)

            reactions_all.append(reactions)

        except:
            reactions_all.append(pd.DataFrame(columns=[sample_name]))

        try:
            secretions = pd.read_table(
                os.path.join(sample_dir, 'secretions.txt'),
                index_col=0)

            secretions_all.append(secretions)

        except:
            secretions_all.append(pd.DataFrame(columns=[sample_name]))

        try:
            uptake = pd.read_table(
                os.path.join(sample_dir, 'uptake.txt'),
                index_col=0)

            uptake_all.append(uptake)

        except:
            uptake_all.append(pd.DataFrame(columns=[sample_name]))

    # Join and output
    reactions_all = pd.concat(reactions_all, axis=1)
    reactions_all.to_csv(
        os.path.join(out_dir, 'reactions.txt'), sep="\t")

    secretions_all = pd.concat(secretions_all, axis=1)
    secretions_all.to_csv(
        os.path.join(out_dir, 'secretions.txt'), sep="\t")

    uptake_all = pd.concat(uptake_all, axis=1)
    uptake_all.to_csv(
        os.path.join(out_dir, 'uptake.txt'), sep="\t")


def load_config(args):
    """
    If a config file is specified, this loads the file
    and applies the arguments in the config file - overwriting
    other arguments specified at the command line

    Really just for batch jobs to make it easier to
    propagate arguments
    """

    if ("config_file" not in args or
            args["config_file"] is None):
        return

    filename = args["config_file"]
    with open(filename) as fin:
        newArgs = json.load(fin)

    args.update(newArgs)
