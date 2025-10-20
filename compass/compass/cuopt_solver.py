from dataclasses import dataclass, field
from build.lib.compass.globals import EXCHANGE_LIMIT
from cuopt.linear_programming.problem import Problem, CONTINUOUS, MAXIMIZE, MINIMIZE
from cuopt.linear_programming.solver_settings import SolverSettings

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

class cuOptSolver:
    """Concrete implementation of LinearProgram using cuOpt solver."""

    def __init__(self, model):
        """Initialize the cuOpt linear program."""
        self.model = model
        self.settings = SolverSettings()
        # Time limit is in seconds.
        self.settings.set_parameter('time_limit', 300)
        self.settings.set_parameter('log_to_console', False)
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
        self.base_problem = None
        # self.solve_problem(LinearProgramDelta(objective=dict(), sense="max"))

    def solve_problem(self, delta: LinearProgramDelta) -> Problem:
        """Creates a linear program based on the model and the provided delta.
        Consult self.problem after calling to get solution status and objective value."""
        problem = Problem("Flux Balance Analysis LP")

        # See if we can apply the delta to the base problem
        reuse_base = len(delta.added_secretion) == 0 and len(delta.added_uptake) == 0 and self.base_problem is not None
        if reuse_base:
            self.problem = self.base_problem

        else:
        variables = {}
        for id, reaction in self.model.reactions.items():
            #print(f"Handling reaction {reaction.id}, {reaction.name}")
            lb = reaction.lower_bound
            ub = reaction.upper_bound
            if id in delta.blocked_reactions:
                ub = reaction.lower_bound
            var = problem.addVariable(lb=lb,
                                             ub=ub,
                                             vtype=CONTINUOUS,
                                             name=reaction.id)
                                             
            variables[id] = var

        # TODO: Double check EXCHANGE_LIMIT vs maximum_flux asymmetry
        # Probably due to limited physical uptake rates vs arbitrary secretion
        for (met_id, rxn_id) in delta.added_secretion.items():
            #print(f"Adding secretion reaction {rxn_id} for metabolite {met_id}")
            secretion_var = self.problem.addVariable(lb=0.0, 
                                            ub=self.model.maximum_flux,
                                            vtype=CONTINUOUS,
                                            name=rxn_id)
            variables[rxn_id] = secretion_var
        for (met_id, rxn_id) in delta.added_uptake.items():
            #print(f"Adding uptake reaction {rxn_id} for metabolite {met_id}")
            uptake_var = self.problem.addVariable(
                                            lb=0.0, 
                                            ub=EXCHANGE_LIMIT,
                                            vtype=CONTINUOUS,
                                            name=rxn_id)
            variables[rxn_id] = uptake_var

        self.variables = variables

        s_mat = self.model.getSMAT()
        constraints = {}
        for metab, rx in s_mat.items():
            #print(f"Handling {metab}: {rx};")
            # Add stoichiometry for added secretion/uptake reactions
            if metab in delta.added_secretion:
                rx.append( (delta.added_secretion[metab], -1.0) )
            if metab in delta.added_uptake:
                rx.append( (delta.added_uptake[metab], 1.0) )
            if len(rx) == 0:
                continue
            
            # rx is a list of (reaction_index, stoichiometric_coeff) tuples
            expr = sum([coeff * variables[id] for id, coeff in rx])

            self.problem.addConstraint(expr == 0, name=metab)
            constraints[metab] = expr

        objective_expr = sum([coeff * variables[id] for id, coeff in delta.objective.items()])
        if delta.sense == "max":
            self.problem.setObjective(objective_expr, sense=MAXIMIZE)
        else:
            self.problem.setObjective(objective_expr, sense=MINIMIZE)
        self.constraints = constraints

    def maximize_reactions(self, reactions: list) -> dict:
        """Maximize the flux through the given reactions using cuOpt."""
        results = {}
        for rxn in reactions:
            blocked_reactions = set()
            objective = {rxn.id: 1.0}
            sense = "max"

            print(f"Maximizing {rxn.id}")

            # Zero out reverse reaction, otherwise extra flux is caused by increasing reverse
            rev_rxn = rxn.reverse_reaction
            if rev_rxn is not None:
                blocked_reactions.add(rev_rxn.id)

            self.solve_problem(LinearProgramDelta(
                blocked_reactions=blocked_reactions,
                objective=objective,
                sense=sense,
            ))
            self.problem.solve(self.settings)

            if self.problem.Status.name == "Optimal":
                print(f"Optimal solution found in {self.problem.SolveTime:.2f} seconds")
                print(f"Objective value = {self.problem.ObjValue}")
                results[rxn.id] = self.problem.ObjValue
            else:
                print(f"Solver ended with status {self.problem.Status.name}")

        return results

    def maximize_metabolites(self, model, metabolites: list) -> dict:
        """Maximize the production of the given metabolites using cuOpt."""
        results = {}
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

            # Note cuopt python does not expose a way to directly get nonzero variables in constraints
            # Though you could iterate over all variables to check
            rxn_ids = self.metab_nonzero[met.id]
            for rxn in [model.reactions[rxn_id] for rxn_id in rxn_ids]:
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

            if secretion_rxn is None:
                added_secretion = True
                secretion_rxn = met.id + "_SECRETION"

                secretion_var = self.problem.addVariable(lb=0.0, 
                                            ub=model.maximum_flux,
                                            vtype=CONTINUOUS,
                                            name=secretion_rxn)






            

            
            




