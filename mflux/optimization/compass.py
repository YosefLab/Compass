"""
Run the procedure for COMPASS
"""
from __future__ import print_function, division, absolute_import
import pulp
from pulp import LpProblem, LpMaximize, LpMinimize, LpAffineExpression
from pulp.solvers import PYGLPK, GLPK_CMD

from . import utils

try:
    from pulp.solvers import CPLEX_PY
except (ImportError):
    print("Error Loading CPLEX")


def run_compass_preprocess(model, debug=False):
    # type: (pandas.Series, mflux.models.MetabolicModel)
    """
    Top-level function for iMat analysis

    Run's COMPASS' preprocessing procedure on a model.
    Results in metabolite and reaction scores.
    """

    # Add reaction fluxes as variables to the problem
    reactions, lbound, ubound, rvars = utils.get_reactions(
        model, exchange_limit=100)

    # Create Constraints
    # Add stoichiometry constraints
    s_constraints = utils.get_connectivity_constraints(model, rvars)

    constraints = s_constraints

    # Create the Problem first
    # Easier to modify existing problem and re-solve
    problem = LpProblem("iMat", LpMaximize)
    for constraint in constraints:
        problem.addConstraint(constraint)

    # For each metabolite, maximize and minimize its prouduction
    # Each metabolite is a constraint.  For each constraint, run the problem
    # Without that constraint, using the constraint as the objective fxn
    m_max = {}
    m_min = {}
    for cc in constraints:
        metabolite = cc.name

        # Remove the constraint from the problem
        _removeConstraint(problem, cc)
        c_dict = {k: v for k, v in cc.items()}
        obj = LpAffineExpression(c_dict)

        # Minimize
        problem.sense = LpMinimize
        problem.objective = obj
        problem.solve(CPLEX_PY(msg=0))
        value = pulp.value(problem.objective)
        m_min[metabolite] = value

        # Maximize
        problem.sense = LpMaximize
        problem.objective = obj
        problem.solve(CPLEX_PY(msg=0))
        value = pulp.value(problem.objective)
        m_max[metabolite] = value

        # Add the constraint back in
        problem.addConstraint(cc)

    # For each reaction, maximize and minimize its flux
    r_max = {}
    r_min = {}
    for rr in reactions:

        problem.objective = {rr: 1}

        # Minimize
        problem.sense = LpMinimize
        problem.solve(CPLEX_PY(msg=0))
        value = pulp.value(problem.objective)
        r_min[metabolite] = value

        # Maximize
        problem.sense = LpMaximize
        problem.solve(CPLEX_PY(msg=0))
        value = pulp.value(problem.objective)
        r_max[metabolite] = value

    return m_min, m_max, r_min, r_max


def _removeConstraint(problem, constraint):
    """
    Removes a constraint from a PuLP problem
    """

    del problem.constraints[constraint.name]
