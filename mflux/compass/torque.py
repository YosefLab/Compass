import os
import pandas as pd
import subprocess as sp

from ..globals import RESOURCE_DIR

TEMPLATE_DIR = os.path.join(RESOURCE_DIR, 'Queue Templates')


def submitCompassTorque(data, model, media, temp_dir, lambda_,
                        output_dir, queue):
    """
    Submits each column of the expression matrix as a distinct job to
    a Torque queueing system

    Parameters
    ==========
    data : str
       Full path to data file

    model : str
        Name of metabolic model to use

    media : str or None
        Name of media to use

    temp_dir : str
        Directory - where to look for sample results.

    out_dir : str
        Where to store aggregated results.  Is created if it doesn't exist.

    lambda_ : float
        Degree of smoothing for single cells.  Valid range from 0 to 1.

    queue : str
        Which queue (name) to submit to
    """

    # Create an array job for single samples
    if media is None:
        media = 'None'

    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)

    # Get the number of samples for array indices
    expression = pd.read_table(data, index_col=0)
    n_samples = len(expression.columns)

    script_args = [data, model, media, lambda_]

    singleSampleScript = os.path.join(TEMPLATE_DIR, "CompassSingleSample.sh")

    command_args = ['qsub', singleSampleScript, '-N', 'COMPASS',
                    '-e', 'localhost:/dev/null',
                    '-o', 'localhost:/dev/null',
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

    # Take the array_job_id and use it to create a collect script
    collectScript = os.path.join(TEMPLATE_DIR, "CompassCollect.sh")
    script_args = [data, model, media, temp_dir]

    command_args = ['qsub', collectScript, '-N', 'COMPASSCollect',
                    '-e', 'localhost:/dev/null',
                    '-o', 'localhost:/dev/null',
                    '-q', queue,
                    '-l', 'nodes=1:ppn=1',
                    '-l', 'walltime=24:00:00',
                    '-l', 'cput=04:00:00',
                    '-V',
                    '-W', 'depend=afterokarray:'+array_job_id,
                    '-d', output_dir,
                    '-F', " ".join(script_args)]

    collect_job_id = sp.check_output(command_args)
    if isinstance(collect_job_id, bytes):
        collect_job_id = collect_job_id.decode()

    print("Compass submitted as array job {} and collect job {}".format(array_job_id, collect_job_id))
    print("Use `qstat -t` to check progress")
