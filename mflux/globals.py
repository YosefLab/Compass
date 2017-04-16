"""
Global variables for use by other modules
"""
import os
import logging
import sys
import multiprocessing

_this_directory = os.path.dirname(os.path.abspath(__file__))

RESOURCE_DIR = os.path.join(_this_directory, "Resources")


# Logging
# Different Logging levels
# Use levels in the logging modules

# How to handle for multiple processes?
# Should look this up
# - Each process has its own logger?
# - Yes - do what Geopy does


logger = logging.getLogger('mflux')
logger.setLevel(logging.DEBUG)

fh = logging.FileHandler("mflux.log", mode='w')
fh.name = 'logfile'

formatter = logging.Formatter(
    '%(asctime)s - %(process)d\n%(message)s\n')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)

logger.addHandler(fh)

# Check if we're the main process
if multiprocessing.current_process().name == 'MainProcess':
    sh = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(message)s\n')
    sh.setFormatter(formatter)
    sh.setLevel(logging.INFO)
    logger.addHandler(sh)
