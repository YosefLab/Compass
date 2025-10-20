from cuopt.linear_programming.problem import Problem, CONTINUOUS, MAXIMIZE
from cuopt.linear_programming.solver_settings import SolverSettings

class cuOptLinearProgram:
    """Concrete implementation of LinearProgram using cuOpt solver."""

    def __init__(self):
        """Initialize the cuOpt linear program."""
        self.problem = Problem("Flux Balance Analysis LP")
        # TODO: Configurable settings
        self.settings = SolverSettings()
        # Time limit is in seconds.
        self.settings.set_parameter('time_limit', 300)
        self.settings.set_parameter('log_to_console', False)
        # Note that there are several tolerance settings, which are generally larger than cplex

    def initialize_problem(self, model):
        """Initialize the cuOpt linear program with the given metabolic model."""

        variables = {}
        for id, reaction in model.reactions.items():
            print(f"Handling reaction {reaction.id}, {reaction.name}")
            var = self.problem.addVariable(lb=reaction.lower_bound,
                                             ub=reaction.upper_bound,
                                             vtype=CONTINUOUS,
                                             name=reaction.id)
            variables[id] = var
        self.variables = variables

        s_mat = model.getSMAT()
        constraints = {}
        for metab, rx in s_mat.items():
            print(f"Handling {metab}: {rx};")
            if len(rx) == 0:
                continue
            
            # rx is a list of (reaction_index, stoichiometric_coeff) tuples
            expr = sum([coeff * variables[id] for id, coeff in rx])
            self.problem.addConstraint(expr == 0, name=metab)

            constraints[metab] = expr
        self.contstraints = constraints


    def maximize_reactions(self, reactions: list) -> dict[int, float]:
        """Maximize the flux through the given reactions using cuOpt."""
        results = {}
        for rxn in reactions:
            print(f"Maximizing {rxn.id}")
            self.problem.setObjective(self.variables[rxn.id], sense = MAXIMIZE)

            # Zero out reverse reaction, otherwise extra flux is caused by increasing reverse
            rev_rxn = rxn.reverse_reaction
            if rev_rxn is not None:
                rev_var = self.variables[rev_rxn.id]  # Fixed typo: varabies -> variables
                old_rev_ub = rev_var.getUpperBound()
                old_rev_lb = rev_var.getLowerBound()
                # Set to 0, or to lower bound, in case lower bound was nonzero
                # This avoids having up = 0 < lb, if the gsmm has unusual bounds
                rev_var.setUpperBound(max(old_rev_lb, 0))

            self.problem.solve(self.settings)

            if self.problem.Status.name == "Optimal":
                print(f"Optimal solution found in {self.problem.SolveTime:.2f} seconds")
                print(f"Objective value = {self.problem.ObjValue}")
                results[rxn.id] = self.problem.ObjValue
            else:
                print(f"Solver ended with status {self.problem.Status.name}")

            if rev_rxn is not None:
                rev_var.setUpperBound(old_rev_ub)

        return results

    def maximize_metabolites(self, metabolites: list[int]) -> dict[int, float]:
        """Maximize the production of the given metabolites using cuOpt."""
        # Implementation specific to cuOpt
        pass
