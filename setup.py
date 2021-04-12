import os
from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy.distutils

# Parse the version string
__version__ = ""
this_directory = os.path.dirname(os.path.abspath(__file__))
version_file = os.path.join(this_directory, "compass", "_version.py")
exec(open(version_file).read())  # Loads version into __version__

# Note that the Cythonization will ONLY effect the Gaussian smoothing which uses the TSNE extension. 
# Extensions
try:
    from Cython.Build import cythonize
    use_cython = True
except ImportError:
    use_cython = False

if use_cython:
    extensions = [
        Extension(
            "compass.compass.extensions.tsne_utils",
            ["compass/compass/extensions/tsne_utils.pyx"],
            include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())
    ]
    extensions = cythonize(extensions)
else:
    extensions = [
        Extension(
            "compass.compass.extensions.tsne_utils",
            ["compass/compass/extensions/tsne_utils.c"],
            include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())
    ]

setup(
    name="compass-sc",
    version=__version__,
    packages=find_packages(),
    ext_modules=extensions,
    include_package_data=True,
    setup_requires = ['numpy>=1.12'], #This has been deprecated, but I like having compatibility with older setuptools versions.
    entry_points={'console_scripts':
                  ['compass = compass.main:entry']},

    install_requires=[
        'pandas>=0.20',
        'tqdm>=4.11',
        'python-libsbml>=5.13',
        'six>=1.10',
        'scikit-learn>=0.19',
        'scipy>=1.0',
        'python-igraph>=0.8',
        'leidenalg'],

        # 'cplex>=12.7.0.0' also required, but installed separately

    author="David DeTomaso",
    author_email="david.detomaso@berkeley.edu",
    description="Metabolic Flux Balance Analysis",
    keywords="",
    url=""
)
