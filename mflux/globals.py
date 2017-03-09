"""
Global variables for use by other modules
"""
import os

_this_directory = os.path.dirname(os.path.abspath(__file__))

RESOURCE_DIR = os.path.join(_this_directory, "Resources")
