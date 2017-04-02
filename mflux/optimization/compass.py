"""
Run the procedure for COMPASS
"""
from __future__ import print_function, division, absolute_import
from tqdm import tqdm

from . import utils
from .. import models

import cplex

__all__ = ['run_compass_preprocess']


def run_compass_preprocess(model):
    # type: (mflux.models.MetabolicModel)
    """
    Cplex-optimized - doesn't use PuLP

    Run's COMPASS' preprocessing procedure on a model.
    Results in metabolite and reaction scores.
    """


    # Add reaction fluxes as variables to the problem
    model.limitUptakeReactions(limit=100)

    # Split fluxes into _pos / _neg
    model.make_unidirectional()

    # Create Constraints
    # Add stoichiometry constraints
    c_lin_expr, c_senses, c_rhs, c_names = (
        utils.get_connectivity_constraints_cplex(model))

    # Index constraints by name for later
    constraints = {}
    for lin_expr, sense, rhs, name in zip(
            c_lin_expr, c_senses, c_rhs, c_names):

        constraints[name] = {
            'lin_expr': lin_expr,
            'sense': sense,
            'rhs': rhs
        }

    # Create the Problem first
    # Easier to modify existing problem and re-solve
    problem = cplex.Cplex()
    problem.set_log_stream(None)  # Suppress output
    problem.set_error_stream(None)  # Suppress errors
    problem.set_warning_stream(None)  # Suppress Warnings
    problem.set_results_stream(None)  # Suppress results to output

    # Add variables
    reactions = list(model.reactions.values())
    problem.variables.add(
        names=[x.id for x in reactions],
        ub=[x.upper_bound for x in reactions],
        lb=[x.lower_bound for x in reactions],)

    # Add constraints
    problem.linear_constraints.add(
        lin_expr=c_lin_expr,
        senses=c_senses,
        rhs=c_rhs,
        names=c_names)

    # For each metabolite, maximize and minimize its prouduction
    # Each metabolite is a constraint.  For each constraint, run the problem
    # Without that constraint, using the constraint as the objective fxn
    m_max = {}
    m_min = {}
    print("COMPASS Preprocess:  Evaluate Metabolite Scores")
    for metabolite in tqdm(constraints):

        # Remove the constraint from the problem
        problem.linear_constraints.delete(metabolite)

        # Set new objective
        sp = constraints[metabolite]['lin_expr']
        ind, val = sp.unpack()
        problem.objective.set_linear(
            [(i, v) for i, v in zip(ind, val)]
        )

        # Minimize
        problem.objective.set_sense(problem.objective.sense.minimize)
        problem.solve()
        value = problem.solution.get_objective_value()
        m_min[metabolite] = value

        # Maximize
        problem.objective.set_sense(problem.objective.sense.maximize)
        problem.solve()
        value = problem.solution.get_objective_value()
        m_max[metabolite] = value

        # Add the constraint back in
        problem.linear_constraints.add(
            lin_expr=[constraints[metabolite]['lin_expr']],
            senses=[constraints[metabolite]['sense']],
            rhs=[constraints[metabolite]['rhs']],
            names=[metabolite],
        )

    # For each reaction, maximize its flux
    r_max = {}
    print("COMPASS Preprocess:  Evaluate Reaction Scores")
    for rr in tqdm(reactions):

        problem.objective.set_linear(
            [(rr, 1)]
        )

        # Maximize
        problem.objective.set_sense(problem.objective.sense.maximize)
        problem.solve()
        value = problem.solution.get_objective_value()
        r_max[metabolite] = value

    return m_min, m_max, r_max
