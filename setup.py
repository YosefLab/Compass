import os
from setuptools import setup, find_packages
from Cython.Build import cythonize
from setuptools.extension import Extension
import numpy.distutils

# Parse the version string
__version__ = ""
this_directory = os.path.dirname(os.path.abspath(__file__))
version_file = os.path.join(this_directory, "mflux", "_version.py")
exec(open(version_file).read())  # Loads version into __version__

# Extensions

extensions = [
    Extension(
        "mflux.compass.extensions.tsne_utils",
        ["mflux/compass/extensions/tsne_utils.pyx"],
        include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())
]

setup(
    name="mflux",
    version=__version__,
    packages=find_packages(),
    ext_modules=cythonize(extensions),
    include_package_data=True,

    entry_points={'console_scripts':
                  ['compass = mflux.compass.main:entry']},

    install_requires=[
        'numpy>=1.12',
        'pandas>=0.20',
        'tqdm>=4.11',
        'python-libsbml>=5.13',
        'six>=1.10',
        'cplex>=12.7.0.0'],

    author="David DeTomaso",
    author_email="david.detomaso@berkeley.edu",
    description="Metabolic Flux Balance Analysis",
    keywords="",
    url=""
)
