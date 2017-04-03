"""
optimization/utils.py

Functions to be used by other optimization routines
"""

from __future__ import print_function, division
from pulp import LpProblem, LpVariable, LpAffineExpression
import cplex

# For pyglpk need to install
# libglpk-dev libglpk36 glpk-utils
# glpk-utils only needed for the CMD version I think

# Also need to install glpk python package
# Might get it here: https://github.com/bradfordboyle/pyglpk, regular one is broken :(
# pip install git+git://github.com/bradfordboyle/pyglpk

# No, the above is not correct
# Need to install python-glpk
# Can do it with sudo apt-get, but then it installs in system python
# Instead, can download it and build it
# Prerequisites include the 'swig' package (sudo apt-get) and
# Instructions here: https://en.wikibooks.org/wiki/GLPK/Python
# Agh, no python 3 support!!!
# Also looks like python-glpk might not work with new versions of glpk


def get_connectivity_constraints(model, rvars):
    """
    Uses the s_mat to define connectivity constraints
    """
    s_constraints = []
    s_mat = model.getSMAT()

    for metab, rx in s_mat.items():
        if len(rx) == 0:
            continue

        exp = LpAffineExpression()
        for reaction, coef in rx:
            exp = exp + rvars[reaction] * coef

        constraint = exp == 0
        constraint.name = metab
        s_constraints.append(constraint)

    return s_constraints


def get_connectivity_constraints_cplex(model):
    """
    Uses the s_mat to define connectivity constraints
    """
    s_mat = model.getSMAT()

    lin_expr = []
    senses = []
    rhs = []
    names = []

    for metab, rx in s_mat.items():
        if len(rx) == 0:
            continue

        ind = [x[0] for x in rx]
        val = [x[1] for x in rx]

        lin_expr.append(cplex.SparsePair(ind=ind, val=val))
        senses.append('E')
        rhs.append(0)
        names.append(metab)

    return lin_expr, senses, rhs, names


def reset_objective(problem):
    """
    Clears all the objective coefficients for the current problem
    by setting all to 0
    """

    names = problem.variables.get_names()
    zeros = [0 for x in names]

    problem.objective.set_linear(zip(names, zeros))
