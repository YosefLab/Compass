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
MODEL_DIR = os.path.join(RESOURCE_DIR, 'Metabolic Models')

# Parameters for COMPASS
BETA = 0.95  # Used to constrain model near optimal point
EXCHANGE_LIMIT = 1.0  # Limit for exchange reactions

# Logging
# Different Logging levels
# Use levels in the logging modules


def init_logger(directory="."):
    logger = logging.getLogger('compass')
    logger.setLevel(logging.DEBUG)
    logger.handlers = []

    # Add file stream to compass.log
    log_file = os.path.join(directory, "compass.log")
    fh = logging.FileHandler(log_file, mode='w')
    fh.name = 'logfile'

    formatter = logging.Formatter("%(levelname)s %(lineno)s %(name)s: %(message)s")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)

    logger.addHandler(fh)

    # Add stream to stdout
    sh = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter("%(levelname)s %(lineno)s %(name)s: %(message)s")
    sh.setFormatter(formatter)
    sh.setLevel(logging.INFO)
    logger.addHandler(sh)
