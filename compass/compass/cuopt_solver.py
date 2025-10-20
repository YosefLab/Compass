from cuopt.linear_programming.problem import Problem, CONTINUOUS, MAXIMIZE
from cuopt.linear_programming.solver_settings import SolverSettings

class cuOptSolver:
    """Concrete implementation of LinearProgram using cuOpt solver."""

    def __init__(self):
        """Initialize the cuOpt linear program."""
        self.settings = SolverSettings()
        # Time limit is in seconds.
        self.settings.set_parameter('time_limit', 300)
        self.settings.set_parameter('log_to_console', False)
        # Note that there are several tolerance settings, which are generally larger than cplex

    # This function must be called first. All following calls should re-establish invariants
    # I.e., reinit problem if you destructively edit the LP
    def initialize_problem(self, model):
        """Initialize the cuOpt linear program with the given metabolic model."""
        self.problem = Problem("Flux Balance Analysis LP")
        # TODO: Configurable settings

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
        metab_nonzero = {}
        for metab, rx in s_mat.items():
            print(f"Handling {metab}: {rx};")
            if len(rx) == 0:
                continue
            
            # rx is a list of (reaction_index, stoichiometric_coeff) tuples
            metab_nonzero[metab] = [id for id, _ in rx]
            expr = sum([coeff * variables[id] for id, coeff in rx])

            self.problem.addConstraint(expr == 0, name=metab)
            constraints[metab] = expr

        self.contstraints = constraints
        # Used for maximize_metabolites
        self.metab_nonzero = metab_nonzero


    def maximize_reactions(self, model, reactions: list) -> dict:
        """Maximize the flux through the given reactions using cuOpt."""
        # Keeping model as parameter for consistency, but this implementation does not need it
        _ = model

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

    def maximize_metabolites(self, model, metabolites: list) -> dict:
        """Maximize the production of the given metabolites using cuOpt."""
        
        used_in_reactions = { c.getConstraintName() for c in self.constraints }

        for met in metabolites:

            if met.id not in used_in_reactions:
                # If the metabolite does not appear in any reactions (this occurs in RECON2)
                # we can skip processing this metabolite.
                continue

            # Find secretion/uptake reactions or create them.
            uptake_rxn = None
            extra_uptake_rxns = []
            secretion_rxn = None
            extra_secretion_rxns = []

            added_uptake = False     # Did we add an uptake reaction?
            added_secretion = False  # "   "   "  "  secretion reaction?

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
                # Then here you need to edit the LP. This does not appear to be something cuopt
                # supports, so it will requre some refactoring.





            

            
            




