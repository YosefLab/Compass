import os
import pandas as pd
import subprocess as sp
import json
import logging

from ..globals import RESOURCE_DIR

TEMPLATE_DIR = os.path.join(RESOURCE_DIR, 'Queue Templates')


def submitCompassTorque(args, temp_dir, output_dir, queue):
    """
    Submits each column of the expression matrix as a distinct job to
    a Torque queueing system

    Parameters
    ==========
    args : dict
        Arguments for the COMPASS call

    temp_dir : str
        Directory - where to look for sample results.

    out_dir : str
        Where to store aggregated results.  Is created if it doesn't exist.

    queue : str
        Which queue (name) to submit to
    """

    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)

    # Get the number of samples for array indices
    data = args['data']
    expression = pd.read_table(data, index_col=0)
    n_samples = len(expression.columns)

    config_file = os.path.join(temp_dir, 'config.json')
    script_args = [config_file]

    # Save the arguments to the config file
    exclude_args = {'temp_dir', 'output_dir', 'torque_queue',
                    'num_processes'}
    newArgs = {}
    for arg, val in args.items():
        if val is None: continue
        if arg in exclude_args: continue
        newArgs[arg] = val

    newArgs['temp_dir'] = '.'

    with open(config_file, 'w') as fp:
        json.dump(newArgs, fp)

    # Submit the single sample array job
    singleSampleScript = os.path.join(TEMPLATE_DIR, "CompassSingleSample.sh")

    command_args = ['qsub', singleSampleScript, '-N', 'COMPASS',
                    #'-e', os.path.join(temp_dir, 'temperr'),
                    #'-o', os.path.join(temp_dir, 'tempout'),
                    '-k', 'oe',
                    '-q', queue,
                    '-l', 'nodes=1:ppn=1',
                    '-l', 'walltime=24:00:00',
                    '-l', 'cput=04:00:00',
                    '-V',
                    '-t', '0-'+str(n_samples-1),
                    '-d', temp_dir,
                    '-F', " ".join(script_args)]

    array_job_id = sp.check_output(command_args)
    if isinstance(array_job_id, bytes):
        array_job_id = array_job_id.decode()

    config_file = os.path.join(temp_dir, 'configCollect.json')
    script_args = [config_file]

    # Save the arguments to the config file
    exclude_args = {'output_dir', 'torque_queue',
                    'num_processes', 'collect'}
    newArgs = {}
    for arg, val in args.items():
        if val is None: continue
        if arg in exclude_args: continue
        newArgs[arg] = val

    newArgs['output_dir'] = '.'

    with open(config_file, 'w') as fp:
        json.dump(newArgs, fp)

    # Take the array_job_id and use it to create a collect script
    collectScript = os.path.join(TEMPLATE_DIR, "CompassCollect.sh")

    command_args = ['qsub', collectScript, '-N', 'COMPASSCollect',
                    '-e', 'localhost:/dev/null',
                    '-o', 'localhost:/dev/null',
                    '-q', queue,
                    '-l', 'nodes=1:ppn=1',
                    '-l', 'walltime=24:00:00',
                    '-l', 'cput=04:00:00',
                    '-V',
                    '-W', 'depend=afteranyarray:'+array_job_id,
                    '-d', output_dir,
                    '-F', " ".join(script_args)]

    collect_job_id = sp.check_output(command_args)
    if isinstance(collect_job_id, bytes):
        collect_job_id = collect_job_id.decode()

    logger = logging.getLogger('mflux')
    logger.info("Compass submitted as array job {} and collect job {}"
                 .format(array_job_id, collect_job_id))
    print("Use `qstat -t` to check progress")
