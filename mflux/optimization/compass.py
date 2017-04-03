"""
Run the procedure for COMPASS
"""
from __future__ import print_function, division, absolute_import
import os
import json
import pandas as pd
from tqdm import tqdm

from . import utils

import cplex

__all__ = ['run_compass_preprocess']

from ..globals import RESOURCE_DIR
PREPROCESS_CACHE_DIR = os.path.join(RESOURCE_DIR, 'COMPASS')
if not os.path.isdir(PREPROCESS_CACHE_DIR):
    os.mkdir(PREPROCESS_CACHE_DIR)

BETA = 0.95  # Used to constrain model near optimal point


def run_compass(model, expression):
    # type: (mflux.models.MetabolicModel, pandas.DataFrame)
    """
    Runs COMPASS on many samples
    """
    EXCHANGE_LIMIT = 100

    # Limit exchange reactions
    model.limitUptakeReactions(limit=EXCHANGE_LIMIT)

    # Split fluxes into _pos / _neg
    model.make_unidirectional()

    # Build model into cplex problem
    problem = build_cplex_problem(model)

    # Run/load pre-processes

    m_uptake, m_secrete, r_max = run_compass_preprocess(model, problem)

    # For each cell/sample
    # Eval expression scores
    # Eval reaction penalties
    # Eval Metabolite secretion/uptake penalties

    all_reaction_scores = {}
    all_uptake_scores = {}
    all_secretion_scores = {}

    for sample in expression.columns:

        expression_data = expression[sample]

        reaction_penalties = eval_reaction_penalties(model, expression_data)

        reaction_scores = compass_reactions(
            model, problem, r_max, reaction_penalties)

        uptake_scores, secretion_scores = compass_exchange(
            model, problem, m_uptake, m_secrete, reaction_penalties)

        all_reaction_scores[sample] = pd.Series(reaction_scores)
        all_uptake_scores[sample] = pd.Series(uptake_scores)
        all_secretion_scores[sample] = pd.Series(secretion_scores)

    reaction_table = pd.DataFrame(all_reaction_scores)
    uptake_table = pd.DataFrame(all_uptake_scores)
    secretion_table = pd.DataFrame(all_secretion_scores)

    return reaction_table, uptake_table, secretion_table


def eval_reaction_penalties(model, expression_data):
    # type: (mflux.models.MetabolicModel, pandas.Series)
    """
    Determines reaction penalties, on the given model, for
    the given expression data
    """

    reaction_expression = model.getReactionExpression(expression_data)
    reaction_expression[reaction_expression < 0] = 0

    reaction_penalties = 1 / (1 + reaction_expression)

    reaction_penalties = reaction_penalties.dropna()

    return reaction_penalties


def compass_exchange(model, problem, m_uptake, m_secrete, reaction_penalties):
    """
    Iterates through metabolites, finding each's max
    uptake and secretion potentials.

    Holds each near its max uptake/secretion while minimizing
    penalty

    Returns the optimal penalty for uptake and secretion

    Returns
    -------
    uptake_scores: dict
        key: species_id
        value: minimum penalty achieved

    secretion_scores: dict
        key: species_id
        value: minimum penalty achieved
    """

    secretion_scores = {}
    uptake_scores = {}
    for metabolite in model.species.values():

        met_id = metabolite.id

        # Rectify exchange reactions
        # Either find existing pos and neg exchange reactions
        # Or create new ones

        uptake_rxn = None
        extra_uptake_rxns = []
        secretion_rxn = None
        extra_secretion_rxns = []

        # Metabolites represented by a constraint: get associated reactions
        sp = problem.linear_constraints.get_rows(met_id)
        rxn_ids = problem.variables.get_names(sp.ind)
        reactions = [model.reaction[x] for x in rxn_ids]

        for reaction in reactions:
            if reaction.is_exchange and met_id in reaction.products:
                if uptake_rxn is None:
                    uptake_rxn = reaction.id
                else:
                    extra_uptake_rxns.append(reaction.id)

            elif reaction.is_exchange and met_id in reaction.reactants:
                if secretion_rxn is None:
                    secretion_rxn = reaction.id
                else:
                    extra_secretion_rxns.append(reaction.id)

        if uptake_rxn is None:
            uptake_rxn = "ADDED_UPTAKE"

            # Add uptake reaction to the problem as a variable
            problem.variables.add(
                names=[uptake_rxn],
                ub=[1000.0],
                lb=[0.0],)

            # Add it to the metabolites constraint
            rxn_index = problem.variables.get_indices(uptake_rxn)
            sp.ind.append(rxn_index)
            sp.val.append(1.0)

        if secretion_rxn is None:
            secretion_rxn = "ADDED_SECRETION"

            # Add secretion reaction to the problem as a variable
            problem.variables.add(
                names=[secretion_rxn],
                ub=[1000.0],
                lb=[0.0],)

            # Add it to the metabolites constraint
            rxn_index = problem.variables.get_indices(secretion_rxn)
            sp.ind.append(rxn_index)
            sp.val.append(-1.0)

        all_uptake = [uptake_rxn] + extra_uptake_rxns
        all_secretion = [secretion_rxn + extra_secretion_rxns]

        # -----------------
        # Optimal Secretion
        # -----------------

        # Close all uptake
        old_uptake_upper = {}
        for rxn_id in all_uptake:
            old_ub = problem.variables.get_upper_bounds(rxn_id)
            old_uptake_upper[rxn_id] = old_ub
            problem.variables.set_upper_bounds(0.0)

        # Close extra secretion
        old_secretion_upper = {}
        for rxn_id in extra_secretion_rxns:
            old_ub = problem.variables.get_upper_bounds(rxn_id)
            old_secretion_upper[rxn_id] = old_ub
            problem.variables.set_upper_bounds(rxn_id, 0.0)

        # # Maximize secretion
        # problem.objective.set_linear(
        #     [(secretion_rxn, 1)]
        # )

        # # Maximize
        # problem.objective.set_sense(problem.objective.sense.maximize)
        # problem.solve()
        # secretion_max = problem.solution.get_objective_value()

        secretion_max = m_secrete[metabolite.id]

        # Set contraint of max secretion with BETA*max
        problem.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=[secretion_rxn], val=[1.0])],
            senses=['R'],
            rhs=[BETA * secretion_max],
            names=['SECRETION_OPT'])

        # Find minimimum penalty
        problem.objective.set_linear(
            list(reaction_penalties.iteritems())
        )

        problem.objective.set_sense(problem.objective.sense.minimize)
        problem.solve()
        value = problem.solution.get_objective_value()
        secretion_scores[met_id] = value

        # Clear Secretion constraint
        problem.linear_constraints.delete('SECRETION_OPT')

        # Restore all uptake
        for rxn_id, old_ub in old_uptake_upper.items():
            problem.variables.set_upper_bounds(rxn_id, old_ub)

        # Restore extra secretion
        for rxn_id, old_ub in old_secretion_upper.items():
            problem.variables.set_upper_bounds(rxn_id, old_ub)

        # -----------------
        # Optimal Uptake
        # -----------------

        # Close extra uptake
        old_uptake_upper = {}
        for rxn_id in extra_uptake_rxns:
            old_ub = problem.variables.get_upper_bounds(rxn_id)
            old_uptake_upper[rxn_id] = old_ub
            problem.variables.set_upper_bounds(0.0)

        # Close all secretion
        old_secretion_upper = {}
        for rxn_id in all_secretion:
            old_ub = problem.variables.get_upper_bounds(rxn_id)
            old_secretion_upper[rxn_id] = old_ub
            problem.variables.set_upper_bounds(rxn_id, 0.0)

        # # Maximize uptake
        # problem.objective.set_linear(
        #     [(uptake_rxn, 1)]
        # )

        # # Maximize
        # problem.objective.set_sense(problem.objective.sense.maximize)
        # problem.solve()
        # uptake_max = problem.solution.get_objective_value()

        uptake_max = m_uptake[metabolite.id]

        # Set contraint of max uptake with BETA*max
        problem.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=[uptake_rxn], val=[1.0])],
            senses=['R'],
            rhs=[BETA * uptake_max],
            names=['UPTAKE_OPT'])

        # Find minimimum penalty
        problem.objective.set_linear(
            list(reaction_penalties.iteritems())
        )

        problem.objective.set_sense(problem.objective.sense.minimize)
        problem.solve()
        value = problem.solution.get_objective_value()
        uptake_scores[met_id] = value

        # Clear Secretion constraint
        problem.linear_constraints.delete('UPTAKE_OPT')

        # Restore extra uptake
        for rxn_id, old_ub in old_uptake_upper.items():
            problem.variables.set_upper_bounds(rxn_id, old_ub)

        # Restore all secretion
        for rxn_id, old_ub in old_secretion_upper.items():
            problem.variables.set_upper_bounds(rxn_id, old_ub)

        # Remove added uptake and secretion reactions
        if uptake_rxn == "ADDED_UPTAKE":
            problem.variables.delete(uptake_rxn)

        if secretion_rxn == "ADDED_SECRETION":
            problem.variables.delete(secretion_rxn)

    return uptake_scores, secretion_scores


def compass_reactions(model, problem, r_max, reaction_penalties):
    """
    Iterates through reactions, holding each near
    its max value while minimizing penalty.

    Minimum overall penalty returned for each reaction

    Returns
    -------
    reaction_scores: dict
        key: reaction id
        value: minimum penalty achieved
    """

    # Create overall penalty as the objective function to minimize
    problem.objective.set_linear(
        list(reaction_penalties.iteritems())
    )

    # Iterate through Reactions
    reaction_scores = {}

    for reaction in model.reactions.values():

        # If Reaction can't carry flux anyways, just continue
        if r_max[reaction.id] == 0:
            reaction_scores[reaction.id] = 0
            continue

        # Get partner reaction or None if it's been pruned already
        if reaction.id.endswith('_pos'):
            partner_id = reaction.id.rsplit("_pos", 1)[0] + "_neg"
            if partner_id in model.reactions:
                partner_reaction = model.reactions[partner_id]
            else:
                partner_reaction = None

        elif reaction.id.endswith('_neg'):
            partner_id = reaction.id.rsplit("_neg", 1)[0] + "_pos"
            if partner_id in model.reactions:
                partner_reaction = model.reactions[partner_id]
            else:
                partner_reaction = None

        else:
            raise Exception("Reaction missing _pos or _neg suffix")

        # Set partner reaction upper-limit to 0 in problem
        # Store old limit for later to restore
        if partner_reaction is not None:
            old_partner_ub = problem.variables.get_upper_bounds(partner_id)
            problem.variables.set_upper_bounds(partner_id, 0.0)

        problem.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=[reaction.id], val=[1.0])],
            senses=['R'],
            rhs=[BETA * r_max[reaction.id]],
            names=['REACTION_OPT'])

        # Minimize Penalty
        problem.objective.set_sense(problem.objective.sense.minimize)
        problem.solve()
        value = problem.solution.get_objective_value()
        reaction_scores[reaction.id] = value

        # Remove Constraint
        problem.linear_constraints.delete('REACTION_OPT')

        # Restore limit of partner reaction to old state
        if partner_reaction is not None:
            problem.variables.set_upper_bounds(partner_id, old_partner_ub)

    return reaction_scores


def build_cplex_problem(model):
    # type: (mflux.models.MetabolicModel)
    """
    Builds and returns a cplex problem representing our metabolic model

    Limits exchange reactions and makes all reactions unidirectional
    by splitting into two components
    """

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

    # Add stoichiometry constraints
    c_lin_expr, c_senses, c_rhs, c_names = (
        utils.get_connectivity_constraints_cplex(model))

    problem.linear_constraints.add(
        lin_expr=c_lin_expr,
        senses=c_senses,
        rhs=c_rhs,
        names=c_names)


    return problem


def preprocess_metabolites(model, problem):
    # type: (mflux.models.MetabolicModel, cplex.Cplex)
    """
    Preprocesses the metabolite limits for secretion and uptake
    """

    # For each metabolite, maximize and minimize its prouduction
    # Each metabolite is a constraint.  For each constraint, run the problem
    # Without that constraint, using the constraint as the objective fxn

    m_max = {} # corresponds to secretion
    m_min = {} # corresponds to uptake
    print("COMPASS Preprocess:  Evaluate Metabolite Scores")
    for metabolite in tqdm(model.species.keys()):

        # Save the constraint to re-add later
        metabolite_constraint = {
            "lin_expr": problem.linear_constraints.get_rows(metabolite),
            "rhs": problem.linear_constraints.get_rhs(metabolite),
            "sense": problem.linear_constraints.get_senses(metabolite),
            "name": metabolite
        }

        # Remove the constraint from the problem
        problem.linear_constraints.delete(metabolite)

        # Set new objective
        sp = metabolite_constraint['lin_expr']
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
            lin_expr=[metabolite_constraint['lin_expr']],
            senses=[metabolite_constraint['sense']],
            rhs=[metabolite_constraint['rhs']],
            names=[metabolite_constraint['name']],
        )

    return m_min, m_max


def preprocess_reactions(model, problem):
    # type: (mflux.models.MetabolicModel, cplex.Cplex)
    """
    Preprocesses the reaction limits
    """

    # For each reaction, maximize its flux
    r_max = {}
    print("COMPASS Preprocess:  Evaluate Reaction Scores")
    for rr in tqdm(model.reactions.keys()):

        problem.objective.set_linear(
            [(rr, 1)]
        )

        # Maximize
        problem.objective.set_sense(problem.objective.sense.maximize)
        problem.solve()
        value = problem.solution.get_objective_value()
        r_max[rr] = value

    return r_max


def run_compass_preprocess(model, problem, use_cache=True):
    # type: (mflux.models.MetabolicModel)
    """
    Cplex-optimized - doesn't use PuLP

    Run's COMPASS' preprocessing procedure on a model.
    Results in metabolite and reaction scores.
    """

    # Check if preprocessed results exist already
    # Load if they do
    cache_file = os.path.join(PREPROCESS_CACHE_DIR, model.name + ".preprocess")

    if os.path.exists(cache_file):
        with open(cache_file) as fin:
            out = json.load(fin)

        m_uptake = out['m_uptake']
        m_secrete = out['m_secrete']
        r_max = out['r_max']

    # Otherwise, run and store results
    m_uptake, m_secrete = preprocess_metabolites(model, problem)

    r_max = preprocess_reactions(model, problem)

    cache_data = {
        'm_uptake': m_uptake,
        'm_secrete': m_secrete,
        'r_max': r_max
    }

    with open(cache_file, 'w') as fout:
        json.dump(cache_data, fout, indent=1)

    return m_uptake, m_secrete, r_max
