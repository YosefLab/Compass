from collections import defaultdict
import logging
import numpy as np
import os
import pandas as pd
import sys
import tempfile
from typing import List, Tuple
from unittest.mock import patch

from turbo_mc.iterative_models.matrix_oracle import MatrixOracle
from turbo_mc.models.exclude_constant_columns_model_wrapper import ExcludeConstantColumnsModelWrapper
from turbo_mc.models.column_normalizer_model_wrapper import ColumnNormalizerModelWrapper
from turbo_mc.models.matrix_completion_fast_als import MatrixCompletionFastALS
from turbo_mc.models.cv_matrix_completion_model import TrainValSplitCVMatrixCompletionModel
from turbo_mc.iterative_models.iterative_matrix_completion_model import IterativeMCMWithGuaranteedSpearmanR2

import compass


logger = logging.getLogger("compass")


class CompassResourceManager():
    r"""
    Hackily exposes some new interfaces in Compass, e.g. determining which cells
    Compass is being run on, the reaction ids, etc.
    Thus, this class knows the details of where Compass output data files are located,
    how to parse them, etc.
    """
    def __init__(self, compass_args: List[str]):
        r"""
        :param compass_args: Arguments used to call Compass (as obtained from sys.argv)
            E.g. ["compass", "--media", "media1", "--data", "compass/Resources/Test Data/rsem_tpmTable_full.txt"]
        """
        self.compass_args = compass_args[:]  # Dependency injection
        with patch.object(sys, 'argv', self.compass_args):
            self.compass_parsed_args = compass.main.parseArgs()

    def get_cell_names(self) -> List[str]:
        r"""
        Returns the list of cell names on which Compass is being run.
        """
        logger.info("CompassResourceManager getting cell names ...")
        data = pd.read_csv(self.compass_parsed_args['data'], sep='\t', nrows=0, index_col=0)
        if "--single-sample" in self.compass_args:
            sample_number = int(self.compass_parsed_args['single_sample'])
            return [data.columns[sample_number]]
        else:
            return list(data.columns)

    def get_reaction_ids(self) -> List[str]:
        r"""
        Returns the list of reaction ids on which Compass is being run.
        """
        logger.info("CompassResourceManager getting reactions ids ...")
        # Run Compass on one cell, patching CPLEX to do nothing.
        with tempfile.TemporaryDirectory() as output_dir_name:
            logger.info(
                f"CompassResourceManager created temporary output directory "
                f"'{output_dir_name}' for calling Compass ...")
            compass_args = self.compass_args[:]
            set_argument(compass_args, "--single-sample", "0")
            set_argument(compass_args, "--output-dir", output_dir_name)

            def cplex_solve_stub(self):
                class DummySolution():
                    def get_objective_value(self):
                        return 0
                self.solution = DummySolution()
            with patch.object(sys, 'argv', compass_args),\
                    patch.object(compass.compass.algorithm.cplex.Cplex, 'solve', new=cplex_solve_stub):
                compass.main.entry()
            reaction_scores = CompassResourceManager(compass_args).get_reaction_scores()
            reaction_ids = list(reaction_scores.index)
            return reaction_ids

    def get_reaction_scores(self) -> pd.DataFrame:
        r"""
        Returns the reaction score DataFrame
        """
        logger.info("CompassResourceManager getting reactions scores ...")
        if "--single-sample" in self.compass_args:
            reaction_scores_path = os.path.join(self.compass_parsed_args['temp_dir'], 'reactions.txt')
        else:
            reaction_scores_path = os.path.join(self.compass_parsed_args['output_dir'], 'reactions.tsv')
        logger.info(f"Trying to read reaction scores from '{reaction_scores_path}' ...")
        reaction_scores = pd.read_csv(reaction_scores_path, sep="\t", index_col=0)
        return reaction_scores

    def get_temporary_directory(self) -> str:
        r"""
        Returns the directory where the temporary data is stored.
        """
        return self.compass_parsed_args['temp_dir']


class CompassOracle(MatrixOracle):
    def __init__(self, compass_args: List[str]):
        r"""
        :param compass_args: Arguments with which Compass is being run.
        """
        logger.info("Initializing CompassOracle ...")
        # Let's first figure out on what cells and reactions Compass is being run on.
        self.compass_args = compass_args[:]
        compass_resource_manager = CompassResourceManager(compass_args=compass_args)
        self.cell_names = np.array(compass_resource_manager.get_cell_names())
        self.reaction_ids = np.array(compass_resource_manager.get_reaction_ids())
        cache_dir = os.path.join(
            compass_resource_manager.get_temporary_directory(),
            "turbo_compass_cache")
        # Write out cell names and reaction ids just for logging / debugging.
        self._save_cell_names_and_reaction_ids_for_debugging(cache_dir)
        # Caching
        self.cache_filepath = os.path.join(cache_dir, "turbo_compass_cache.txt")
        # Just make sure that the cache works before even starting!
        self._cache_observations([], [], [], self.cache_filepath)
        # Load the cache
        self._load_cache()

    def _load_cache(self):
        logger.info("CompassOracle trying to load its cache ...")
        rows, cols, vals = [], [], []
        with open(self.cache_filepath, "r") as cache_file:
            for line in cache_file:
                row, col, val = line.strip('\n').split(' ')
                rows.append(int(row)), cols.append(int(col)), vals.append(float(val))
        logger.info(f"CompassOracle loaded {len(vals)} observations from the cache!")
        self._add_observations(rows, cols, vals)

    def _save_cell_names_and_reaction_ids_for_debugging(self, cache_dir: str) -> None:
        r"""
        This is only used for logging / debugging purposes.
        If the cache already exists, will validate its contents.
        """
        cell_names_filepath = os.path.join(cache_dir, "cell_names.txt")
        reaction_ids_filepath = os.path.join(cache_dir, "reaction_ids.txt")

        def validate_cache_and_write_out_list_of_strings(filepath, names: List[str]) -> None:
            r"""
            Writes out the list of strings 'names' to filepath, with one string per line.
            First checks to see if the names are already cached. If they are, checks
            that they are correct.
            """
            # If the cache already exists, check that it has the right content!
            # Or else, something is VERY wrong with the cache and we should raise.
            output_str = '\n'.join(names) + '\n'
            if os.path.exists(filepath):
                with open(filepath) as f:
                    file_contents = f.read()
                    if not file_contents == output_str:
                        raise Exception(f"CompassOracle cache file '{filepath}' seems to be corrupted!")
                    else:
                        # All is well! No need to write anything.
                        logger.info(f"CompassOracle cache file '{filepath}' is healthy!")
                        return
            # Cache does not exist yet, create.
            logger.info(f"Creating CompassOracle cache file '{filepath}' ...")
            if not os.path.exists(os.path.dirname(filepath)):
                os.makedirs(os.path.dirname(filepath))
            with open(filepath, "w") as cache_file:
                cache_file.write(output_str)
        validate_cache_and_write_out_list_of_strings(cell_names_filepath, list(self.cell_names))
        validate_cache_and_write_out_list_of_strings(reaction_ids_filepath, list(self.reaction_ids))

    def shape(self) -> Tuple[int, int]:
        r"""
        Simply returns the shape of the cell X reaction matrix.
        """
        return len(self.cell_names), len(self.reaction_ids)

    def _observe_entries(
        self,
        rows: List[int],
        cols: List[int]
    ) -> np.array:
        # Start by running Compass on the selected cells (rows) and reactions (cols).
        assert("--selected-reactions-for-each-cell" not in self.compass_args)
        compass_args_to_observe_entries = self.compass_args[:]
        with tempfile.NamedTemporaryFile("w") as selected_reactions_file:
            # Create list of selected cells and reactions
            logger.info("CompassOracle populating the file with the selected cells and reactions ...")
            self._populate_selected_reactions_file(rows, cols, selected_reactions_file.name)
            # Now run Compass
            set_argument(
                compass_args_to_observe_entries,
                "--selected-reactions-for-each-cell",
                selected_reactions_file.name)
            with tempfile.TemporaryDirectory() as output_dir_name:
                logger.info(
                    f"CompassOracle created temporary output directory "
                    f"'{output_dir_name}' for running Compass")
                set_argument(compass_args_to_observe_entries, "--output-dir", output_dir_name)
                logger.info("CompassOracle running Compass for the selected subset of the reaction score matrix ...")
                with patch.object(sys, 'argv', compass_args_to_observe_entries):
                    compass.main.entry()
                logger.info(
                    "CompassOracle successfully ran Compass. "
                    "Now will retrieve the computed reaction scores ...")
                # Retrieve reaction scores.
                compass_resource_manager = CompassResourceManager(compass_args_to_observe_entries)
                reaction_scores = compass_resource_manager.get_reaction_scores()
                # Check that rows and cols agree
                assert(list(reaction_scores.columns) == list(self.cell_names))
                assert(list(reaction_scores.index) == list(self.reaction_ids))
                vals = reaction_scores.to_numpy().T[(rows, cols)]
                logger.info("CompassOracle caching the latest Compass results ...")
                self._cache_observations(rows, cols, list(vals), self.cache_filepath)
                return vals

    def _populate_selected_reactions_file(
        self,
        rows: List[int],
        cols: List[int],
        filename: str
    ) -> None:
        r"""
        Populates the file which will be later passed to Compass via the
        --selected-reactions-for-each-cell argument.
        """
        with open(filename, "w") as file:
            row_to_cols_mapping = defaultdict(list)
            for (r, c) in zip(rows, cols):
                row_to_cols_mapping[r].append(c)
            for row, columns in row_to_cols_mapping.items():
                cell_name = self.cell_names[row]
                reaction_ids = self.reaction_ids[columns]
                tokens = [cell_name] + list(reaction_ids)
                line_for_this_cell = ','.join(tokens) + '\n'
                file.write(line_for_this_cell)

    @staticmethod
    def _cache_observations(
            rows: List[int],
            cols: List[int],
            vals: List[float],
            cache_filepath: str):
        r"""
        Adds the newly observed (row, col, val) triples to the cache. The structure of
        the cache is one line per observation, e.g.:
        1 5 1356.7
        37 134 0.0
        123 4567 8901.2
        ...
        """
        assert(len(rows) == len(cols))
        assert(len(rows) == len(vals))
        output_str = ""
        for row, col, val in zip(rows, cols, vals):
            output_str += f"{row} {col} {val}\n"
        if not os.path.exists(os.path.dirname(cache_filepath)):
            os.makedirs(os.path.dirname(cache_filepath))
        with open(cache_filepath, "a+") as cache_file:
            cache_file.write(output_str)


def set_argument(args: List, arg_name: str, val: str) -> None:
    r"""
    Modifies args list in place, setting 'arg_name' argument to 'val'.
    If 'arg_name' does not exist already, it is created.
    """
    if not isinstance(arg_name, str) or not isinstance(val, str):
        raise ValueError(f"{arg_name} and {val} should be strings!")  # pragma: no cover
    for i in range(len(args) - 1):
        if args[i] == arg_name:
            args[i + 1] = str(val)
            return
    args.append(arg_name)
    args.append(val)


def pop_argument(args: List[str], arg: str) -> Tuple[List[str], str]:
    r"""
    Removes the argument from the parameter list and returns the new parameter
    list and the popped parameter's value. If the argument is not present, returns
    the original arguments and the empty string.
    e.g. pop_argument(['compass', '--turbo', '0.95']) == (['compass'], '0.95')
    after the operation.
    """
    if arg not in args:
        return args[:], ""
    for i in range(len(args) - 1):
        if args[i] == arg:
            return args[:i] + args[(i + 2):], args[i + 1]
    assert(False)


def turbo_compass_entry() -> None:
    logger.info("******** Turbo Compass started ********")
    # Extract turbo arguments.
    compass_parsed_args = compass.main.parseArgs()
    requested_cv_spearman_r2, increments, min_pct_meet_sr2_requirement, max_iters =\
        compass_parsed_args['turbo'],\
        compass_parsed_args['turbo_increments'],\
        compass_parsed_args['turbo_min_pct'],\
        compass_parsed_args['turbo_max_iters']
    logger.info(
        f"Running Turbo Compass with arguments:\n"
        f"--turbo={requested_cv_spearman_r2}\n"
        f"--turbo-increments={increments}\n"
        f"--turbo-min-pct={min_pct_meet_sr2_requirement}\n"
        f"--turbo-max-iters={max_iters}")

    # Remove turbo arguments from the command line call, and create the CompassOracle for the new command line.
    compass_args = sys.argv[:]
    compass_args, _ = pop_argument(compass_args, '--turbo')
    compass_args, _ = pop_argument(compass_args, '--turbo-increments')
    compass_args, _ = pop_argument(compass_args, '--turbo-min-pct')
    compass_args, _ = pop_argument(compass_args, '--turbo-max-iters')
    logger.info("Turbo Compass creating the CompassOracle ...")
    compass_matrix_oracle = CompassOracle(compass_args[:])

    model =\
        IterativeMCMWithGuaranteedSpearmanR2(
            cv_model=lambda state:
                TrainValSplitCVMatrixCompletionModel(
                    ExcludeConstantColumnsModelWrapper(
                        ColumnNormalizerModelWrapper(
                            MatrixCompletionFastALS(
                                n_factors=max(1, int(min(state.R, state.C) * state.sampled_density * 0.5)),
                                lam=10.0,
                                n_epochs=100,
                                verbose=False))),
                    train_ratio=0.8,
                    verbose=True,
                    finally_refit_model=None),
            requested_cv_spearman_r2=requested_cv_spearman_r2,
            sampling_density=increments,
            finally_refit_model=lambda state:
                ExcludeConstantColumnsModelWrapper(
                    ColumnNormalizerModelWrapper(
                        MatrixCompletionFastALS(
                            n_factors=int(min(state.R, state.C) * state.sampled_density * 0.5),
                            lam=10.0,
                            n_epochs=300,
                            verbose=False))),
            min_pct_meet_sr2_requirement=min_pct_meet_sr2_requirement,
            verbose=True,
            plot_progress=False,
            max_iterations=max_iters
        )

    # Fit iterative model. This is the core procedure of Turbo Compass (here is where all the magic happens)
    logger.info("Turbo Compass performing iterative matrix completion ...")
    np.random.seed(1)
    model.fit(compass_matrix_oracle)
    X_imputed = model.impute_all().T
    logger.info("Iterative matrix completion was successfull!")

    # Write out imputation
    logger.info("Turbo Compass writing out imputation ...")
    if "--single-sample" in compass_args:
        if not os.path.isdir(compass_parsed_args['temp_dir']):
            os.makedirs(compass_parsed_args['temp_dir'])
        reaction_scores_path = os.path.join(compass_parsed_args['temp_dir'], 'reactions.txt')
    else:
        if not os.path.isdir(compass_parsed_args['output_dir']):
            os.makedirs(compass_parsed_args['output_dir'])
        reaction_scores_path = os.path.join(compass_parsed_args['output_dir'], 'reactions.tsv')
    cell_names = compass_matrix_oracle.cell_names[:]
    reaction_ids = compass_matrix_oracle.reaction_ids[:]
    X_imputed_df = pd.DataFrame(data=X_imputed, index=reaction_ids, columns=cell_names)
    X_imputed_df.to_csv(reaction_scores_path, sep="\t")
    logger.info("Turbo Compass completed successfully!")
    return
