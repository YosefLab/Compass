"""
Global variables for use by other modules
"""
import os
import logging
import sys

_this_directory = os.path.dirname(os.path.abspath(__file__))

GIT_DIR = os.path.abspath(
    os.path.join(_this_directory, '..', '.git')
)

RESOURCE_DIR = os.path.join(_this_directory, "Resources")

NUM_THREADS = 1   # Threads per sample
TEST_MODE = False  # When true, only a small subset of reactions are used

# Logging
# Different Logging levels
# Use levels in the logging modules


def init_logger(directory="."):
    logger = logging.getLogger('mflux')
    logger.setLevel(logging.DEBUG)
    logger.handlers = []

    # Add file stream to mflux.log
    log_file = os.path.join(directory, "mflux.log")
    fh = logging.FileHandler(log_file, mode='w')
    fh.name = 'logfile'

    formatter = logging.Formatter(
        '%(message)s')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)

    logger.addHandler(fh)

    # Add stream to stdout
    sh = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(message)s')
    sh.setFormatter(formatter)
    sh.setLevel(logging.INFO)
    logger.addHandler(sh)
