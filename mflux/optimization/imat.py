"""
Run the procedure for iMAT
"""
from __future__ import print_function, division, absolute_import
import numpy as np
import pulp
import time
from pulp import LpProblem, LpMaximize, LpVariable, LpAffineExpression
from pulp.solvers import PYGLPK, GLPK_CMD

from . import utils

try:
    from pulp.solvers import CPLEX_PY
except (ImportError):
    print("Error Loading CPLEX")


def run_iMat(exp_data, model, debug=False):
    # type: (pandas.Series, mflux.models.MetabolicModel)
    """
    Top-level function for iMat analysis

    Run's iMat on a dataset,  yielding a set of flux predictions
    """
    t_start = time.time()
    # Define iMAT constants
    RH_THRESH = 6
    RL_THRESH = 2
    EPS = 1

    # Add reaction fluxes as variables to the problem
    reactions, lbound, ubound, rvars = utils.get_reactions(model, exchange_limit=100)

    # Score the reactions
    scores = model.getReactionScores(exp_data)

    # Partition the reactions into High and Low
    high_reactions = set()
    low_reactions = set()

    for reaction in scores:

        if np.isnan(scores[reaction]): continue

        if scores[reaction] > RH_THRESH:
            high_reactions.add(reaction)
        elif scores[reaction] < RL_THRESH:
            low_reactions.add(reaction)

    # Add boolean variables to the problem
    # For 'high' reactions (2 var's each)
    rh_yp = {}
    rh_ym = {}
    for reaction in high_reactions:
        rvar_yp = LpVariable(reaction + "_yp", cat=pulp.LpBinary)
        rvar_ym = LpVariable(reaction + "_ym", cat=pulp.LpBinary)
        rh_yp[reaction] = rvar_yp
        rh_ym[reaction] = rvar_ym

    # For 'low' reactions (1 var)
    rl_y = {}
    for reaction in low_reactions:
        rvar_y = LpVariable(reaction + "_y", cat=pulp.LpBinary)
        rl_y[reaction] = rvar_y

    # Create Constraints
    # Add stoichiometry constraints
    s_constraints = utils.get_connectivity_constraints(model, rvars)

    # Eq. 3 constraints (from Paper)
    # For each r in rh (high reactions)
    # v + yi_p(vmin,i - eps) >= vmin,i
    eq_3_constraints = []
    for reaction in high_reactions:
        v = rvars[reaction]
        yp = rh_yp[reaction]
        vmin = lbound[reaction]
        constraint = v + yp * (vmin - EPS) >= vmin
        eq_3_constraints.append(constraint)

    # Eq. 4 constraints (from Paper)
    # For each r in rh (high reactions)
    # v + yi_m(vmax,i + eps) <= vmax,i
    eq_4_constraints = []
    for reaction in high_reactions:
        v = rvars[reaction]
        ym = rh_ym[reaction]
        vmax = ubound[reaction]
        constraint = v + ym * (vmax + EPS) <= vmax
        eq_4_constraints.append(constraint)

    # Eq. 5 constraints (from Paper)
    # For each r in rl (low reactions)
    # v,i >= vmin,i*(1-y,i)
    # v,i <= vmax,i(1-y,i)
    eq_5_constraints = []
    for reaction in low_reactions:
        v = rvars[reaction]
        y = rl_y[reaction]
        vmax = ubound[reaction]
        vmin = lbound[reaction]

        c1 = v >= vmin * (1 - y)
        c2 = v <= vmax * (1 - y)
        eq_5_constraints.append(c1)
        eq_5_constraints.append(c2)

    constraints = s_constraints + eq_3_constraints + eq_4_constraints + eq_5_constraints

    # Create the objective function
    obj = LpAffineExpression()
    for reaction in high_reactions:
        yp = rh_yp[reaction]
        ym = rh_ym[reaction]
        obj = obj + (yp + ym)

    for reaction in low_reactions:
        y = rl_y[reaction]
        obj = obj + y

    # Solve equation
    problem = LpProblem("iMat", LpMaximize)
    for constraint in constraints:
        problem.addConstraint(constraint)

    problem += obj

    t_prep = time.time()
    # problem.solve(GLPK_CMD(msg=0))
    problem.solve(CPLEX_PY(msg=0))
    
    t_solved = time.time()

    # Gather results and return
    result_fluxes = {}
    for reaction in reactions:
        r_val = pulp.value(rvars[reaction])
        result_fluxes[reaction] = r_val
        
    t_output = time.time()

    if debug:
        dt_load = t_prep - t_start
        dt_solve = t_solved-t_prep
        dt_output = t_output-t_solved
        print("Time elapsed: \n  Load:", dt_load, "\n  Solve: ", dt_solve, "\n  Output: ", dt_output)

    return result_fluxes, problem
