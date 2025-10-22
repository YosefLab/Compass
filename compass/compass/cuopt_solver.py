from dataclasses import dataclass, field
import logging
import gc
import cupy

import numpy as np
from build.lib.compass.globals import EXCHANGE_LIMIT
from cuopt.linear_programming.problem import Problem, CONTINUOUS, MAXIMIZE, MINIMIZE
from cuopt.linear_programming.solver_settings import SolverSettings

logger = logging.getLogger(__name__)

@dataclass
class LinearProgramDelta:
    """Track deltas to the GSMM-derived linear program.
    This enables changes that a solver does not support in place, like editing constraints."""
    # Reaction ID -> coefficient
    objective: dict
    # Sense of optimization: max or min
    sense: str
    # Metabolite ID -> reaction ID
    added_secretion: dict = field(default_factory=dict)
    # Metabolite ID -> reaction ID
    added_uptake: dict = field(default_factory=dict)
    # Reaction IDs
    blocked_reactions: set = field(default_factory=set)

@dataclass
class Solution:
    """Solution from a linear program"""
    obj_status: str
    obj_value: np.float64

class cuOptSolver:
    """Concrete implementation of LinearProgram using cuOpt solver."""

    def __init__(self, model):
        """Initialize the cuOpt linear program."""
        self.model = model
        self.settings = SolverSettings()
        # Time limit is in seconds.
        self.settings.set_parameter('time_limit', 300)
        # Note https://github.com/NVIDIA/cuopt/issues/187 means the setting still get sent to stdout
        self.settings.set_parameter('log_to_console', False)
        self.settings.set_parameter('presolve', True)
        # Note that there are several tolerance settings, which are generally larger than cplex
        
        s_mat = self.model.getSMAT()
        metab_nonzero = {}
        for metab, rx in s_mat.items():
            if len(rx) == 0:
                continue
            # rx is a list of (reaction_index, stoichiometric_coeff) tuples
            metab_nonzero[metab] = [id for id, _ in rx]
        self.metab_nonzero = metab_nonzero

        # Create a base problem that can be modified for some deltas
        # For cuopt, we can easily change the upper/lower bounds and objective
        # Changing constraints is much harder, so for any metabolites we reconstruct the whole problem
        self.base_problem = cuOptSolver.create_problem(self.model, LinearProgramDelta(objective=dict(), sense="max"))

    def solve_problem(self, delta: LinearProgramDelta) -> Solution:
        """Solves a linear program based on the model and the provided delta."""

        # See if we can apply the delta to the base problem
        reuse_base = len(delta.added_secretion) == 0 and len(delta.added_uptake) == 0
        restore_upper_bounds = {}
        
        if reuse_base:
            problem = self.base_problem
            for rxn_id in delta.blocked_reactions:
                var = problem.getVariable(rxn_id)
                restore_upper_bounds[rxn_id] = var.getUpperBound()
                var.setUpperBound(var.getLowerBound())
            objective_expr = sum([coeff * problem.getVariable(id) for id, coeff in delta.objective.items()])
            if delta.sense == "max":
                problem.setObjective(objective_expr, sense=MAXIMIZE)
            else:
                problem.setObjective(objective_expr, sense=MINIMIZE)

        else:
            # Create a new problem from scratch
            # TODO: Support editing constraints to reuse problem for metabolites
            problem = cuOptSolver.create_problem(self.model, delta)
 
        try:
            problem.solve(self.settings)
        except Exception as e:
            logger.error(f"cuOPT solve failed!")
            logger.error(f"Error: {type(e).__name__}: {e}")
            raise e

        sol = Solution(obj_status=problem.Status.name, obj_value=problem.ObjValue)

        # Restore problem's constraints
        for (rxn_id, ub) in restore_upper_bounds.items():
            problem.getVariable(rxn_id).setUpperBound(ub)
        
        # Clean up created problem
        # I am optimistic this will limit the observed increasing memory usage
        if not reuse_base:
            del problem
        
        return sol
        
    @staticmethod
    def create_problem(model, delta: LinearProgramDelta):
        """Creates a linear program based on the model and the provided delta."""
        problem = Problem("Flux Balance Analysis LP")
        variables = {}
        for id, reaction in model.reactions.items():
            lb = reaction.lower_bound
            ub = reaction.upper_bound
            if id in delta.blocked_reactions:
                ub = reaction.lower_bound
            var = problem.addVariable(lb=lb, ub=ub, vtype=CONTINUOUS, name=reaction.id)       
            variables[id] = var

        # TODO: Double check EXCHANGE_LIMIT vs maximum_flux asymmetry
        # Probably due to limited physical uptake rates vs arbitrary secretion
        for (met_id, rxn_id) in delta.added_secretion.items():
            secretion_var = problem.addVariable(lb=0.0, 
                                            ub=model.maximum_flux,
                                            vtype=CONTINUOUS,
                                            name=rxn_id)
            variables[rxn_id] = secretion_var
        for (met_id, rxn_id) in delta.added_uptake.items():
            uptake_var = problem.addVariable(
                                            lb=0.0, 
                                            ub=EXCHANGE_LIMIT,
                                            vtype=CONTINUOUS,
                                            name=rxn_id)
            variables[rxn_id] = uptake_var

        s_mat = model.getSMAT()
        constraints = {}
        for metab, rx in s_mat.items():
            # Add stoichiometry for added secretion/uptake reactions
            if metab in delta.added_secretion:
                rx.append( (delta.added_secretion[metab], -1.0) )
            if metab in delta.added_uptake:
                rx.append( (delta.added_uptake[metab], 1.0) )
            if len(rx) == 0:
                continue
            
            # rx is a list of (reaction_index, stoichiometric_coeff) tuples
            expr = sum([coeff * variables[id] for id, coeff in rx])

            problem.addConstraint(expr == 0, name=metab)
            constraints[metab] = expr

        objective_expr = sum([coeff * variables[id] for id, coeff in delta.objective.items()])
        if delta.sense == "max":
            problem.setObjective(objective_expr, sense=MAXIMIZE)
        else:
            problem.setObjective(objective_expr, sense=MINIMIZE)
            
        return problem
        

    def maximize_reactions(self, reactions: list) -> dict:
        """Maximize the flux through the given reactions using cuOpt."""
        results = {}
        solve_count = 0
        
        for rxn in reactions:
            blocked_reactions = set()
            objective = {rxn.id: 1.0}
            sense = "max"

            # Zero out reverse reaction, otherwise extra flux is caused by increasing reverse
            rev_rxn = rxn.reverse_reaction
            if rev_rxn is not None:
                blocked_reactions.add(rev_rxn.id)

            sol = self.solve_problem(LinearProgramDelta(
                blocked_reactions=blocked_reactions,
                objective=objective,
                sense=sense,
            ))

            if sol.obj_status == "Optimal":
                logger.info(f"Reaction {rxn.id}: Objective value = {sol.obj_value};")
                results[rxn.id] = sol.obj_value
            else:
                logger.info(f"Reaction {rxn.id}: Solver ended with status {sol.obj_status}")
            
            # Periodic garbage collection to prevent memory buildup
            solve_count += 1
            if solve_count % 100 == 0:
                logger.debug(f"Completed {solve_count} reactions, triggered GC")
                self.cleanup_memory_and_gpu()

        return results

    def maximize_metabolites(self, metabolites: list) -> dict:
        """Maximize the production of the given metabolites using cuOpt."""
        results = {}
        solve_count = 0
        
        for met in metabolites:

            if met.id not in self.metab_nonzero:
                # If the metabolite does not appear in any reactions (this occurs in RECON2)
                # we can skip processing this metabolite.
                continue

            # Find secretion/uptake reactions or create them.
            uptake_rxn = None
            extra_uptake_rxns = []
            secretion_rxn = None
            extra_secretion_rxns = []
            # For optimal secretion, block uptake and other secretion reactions
            secretion_blocked = set()
            uptake_blocked = set()

            # Note cuopt python does not expose a way to directly get nonzero variables in constraints
            # Though you could iterate over all variables to check
            rxn_ids = self.metab_nonzero[met.id]
            for rxn in [self.model.reactions[rxn_id] for rxn_id in rxn_ids]:
                if rxn.is_exchange and met.id in rxn.products:
                    if uptake_rxn is None:
                        uptake_rxn = rxn.id
                    else:
                        extra_uptake_rxns.append(rxn.id)

                elif rxn.is_exchange and met.id in rxn.reactants:
                    if secretion_rxn is None:
                        secretion_rxn = rxn.id
                    else:
                        extra_secretion_rxns.append(rxn.id)

            added_secretion = {}
            if secretion_rxn is None:
                secretion_rxn = met.id + "_SECRETION"
                logger.debug(f"Adding {secretion_rxn}")
                added_secretion[met.id] = secretion_rxn
            else:
                uptake_blocked.add(secretion_rxn)

            added_uptake = {}
            if uptake_rxn is None:
                uptake_rxn = met.id + "_UPTAKE"
                logger.debug(f"Adding {uptake_rxn}")
                added_uptake[met.id] = uptake_rxn
            else:
                secretion_blocked.add(uptake_rxn)
                
            # Optimal secretion
            # Close all uptake reactions and extra secretion reactions
            blocked_reactions = secretion_blocked | set(extra_uptake_rxns + extra_secretion_rxns)
            objective = { secretion_rxn: 1.0 }
            sense = "max"

            sol = self.solve_problem(LinearProgramDelta(
                objective=objective,
                sense=sense,
                blocked_reactions=blocked_reactions,
                added_secretion=added_secretion,
            ))
            
            if sol.obj_status == "Optimal":
                logger.info(f"Secretion {secretion_rxn}: Objective value = {sol.obj_value}")
                results[secretion_rxn] = sol.obj_value
            else:
                logger.warning(f"Secretion {secretion_rxn}: Solver ended with status {sol.obj_status}")

            # Optimal uptake
            # Close all secretion and extra uptake reactions
            blocked_reactions = uptake_blocked | set(extra_uptake_rxns + extra_secretion_rxns)
            objective = { uptake_rxn: 1.0 }
            sense = "max"

            sol = self.solve_problem(LinearProgramDelta(
                objective=objective,
                sense=sense,
                blocked_reactions=blocked_reactions,
                added_uptake=added_uptake,
            ))

            if sol.obj_status == "Optimal":
                logger.info(f"Uptake {uptake_rxn}: Objective value = {sol.obj_value}")
                results[uptake_rxn] = sol.obj_value
            else:
                logger.warning(f"Uptake {uptake_rxn}: Solver ended with status {sol.obj_status}")
            
            # Periodic garbage collection to prevent memory buildup
            solve_count += 1  #
            if solve_count % 50 == 0:
                logger.debug(f"Completed {solve_count} metabolites, triggered GC")
                self.cleanup_memory_and_gpu()

        return results
    
    def cleanup_memory_and_gpu(self):
        """Force cleanup of GPU memory. Call this periodically during long runs."""
        try:
            # Synchronize all GPU operations
            cupy.cuda.Device(0).synchronize()
            # Clear memory pool
            mempool = cupy.get_default_memory_pool()
            mempool.free_all_blocks()
            logger.debug(f"GPU memory pool cleared")
        except Exception as e:
            logger.warning(f"Error during GPU cleanup: {e}")
        
        # Force Python garbage collection
        gc.collect()
        logger.debug("Python GC triggered")






            

            
            




