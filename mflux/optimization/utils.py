"""
optimization/utils.py

Functions to be used by other optimization routines
"""

from __future__ import print_function, division
from pulp import LpProblem, LpVariable, LpAffineExpression

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
    ex_metabolites = model.getExtracellularMetabolites()

    for metab, rx in s_mat.items():
        if len(rx) == 0:
            continue

        if metab in ex_metabolites:
            continue

        exp = LpAffineExpression()
        for reaction, coef in rx:
            exp = exp + rvars[reaction] * coef

        constraint = exp == 0
        constraint.name = metab
        s_constraints.append(constraint)

    return s_constraints


def get_reactions(model, exchange_limit=None):
    """
    Uses the bounds in the model (global bounds) to define lpVariables

    `exchange_limit` is the limit that is additionally imposed on exchange
    reactions with the extracellular environment

    """
    reactions = model.getReactions()
    lbound, ubound = model.getReactionBounds()

    if exchange_limit is not None:
        lbound, ubound = model.limitUptakeReactions(lbound, ubound, exchange_limit)

    rvars = {}
    for reaction in reactions:
        rvar = LpVariable(reaction, lowBound=lbound[reaction], upBound=ubound[reaction])
        rvars[reaction] = rvar

    return reactions, lbound, ubound, rvars


