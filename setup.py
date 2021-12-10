import os
from setuptools import setup, find_packages

# Parse the version string
__version__ = ""
this_directory = os.path.dirname(os.path.abspath(__file__))
version_file = os.path.join(this_directory, "compass", "_version.py")
exec(open(version_file).read())  # Loads version into __version__

setup(
    name="compass-sc",
    version=__version__,
    packages=find_packages(),
    include_package_data=True,

    entry_points={'console_scripts':
                  ['compass = compass.main:entry']},

    install_requires=[
        'numpy>=1.12',
        'pandas>=0.20',
        'tqdm>=4.11',
        'python-libsbml>=5.13',
        'six>=1.10',
        'scikit-learn>=0.19',
        'scipy>=1.0',
        'python-igraph>=0.9', #Need the version to be high enough fgo
        'leidenalg>=0.8.2',
        'anndata'],

        # 'cplex>=12.7.0.0' also required, but installed separately

    author="David DeTomaso",
    author_email="david.detomaso@berkeley.edu",
    description="Metabolic Flux Balance Analysis",
    keywords="",
    url=""
)
