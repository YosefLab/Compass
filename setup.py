import os
from setuptools import setup, find_packages

# Parse the version string
__version__ = ""
this_directory = os.path.dirname(os.path.abspath(__file__))
version_file = os.path.join(this_directory, "mflux", "_version.py")
exec(open(version_file).read())  # Loads version into __version__

setup(
    name="mflux",
    version=__version__,
    packages=find_packages(),

    install_requires=[],

    author="David DeTomaso",
    author_email="david.detomaso@berkeley.edu",
    description="Metabolic Flux Balance Analysis",
    keywords="",
    url=""
)
