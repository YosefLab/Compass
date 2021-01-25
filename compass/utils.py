"""
utils.py

Functions to be used by other optimization routines
"""

from __future__ import print_function, division
import cplex


def get_steadystate_constraints(model):
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
    problem.objective.set_name('none')
