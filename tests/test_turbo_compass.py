import hashlib
import numpy as np
import os
import pandas as pd
import pytest
import sys
import tempfile
from unittest.mock import patch

import compass.global_state as global_state
import compass.main
from compass.turbo_compass import CompassResourceManager, CompassOracle, set_argument,\
    pop_argument


COMPASS_ARGS = [
    "compass",
    "--media", "media1",
    "--data", "compass/Resources/Test Data/rsem_tpmTable_full.txt"
]


def test_set_argument():
    args = ['compass', '--turbo', '0.95']
    set_argument(args, '--turbo', '1.0')
    assert(args == ['compass', '--turbo', '1.0'])
    args = ['compass']
    set_argument(args, '--turbo', '0.6')
    assert(args == ['compass', '--turbo', '0.6'])


def test_pop_argument():
    args = ['compass', '--turbo', '0.95']
    assert(pop_argument(args, '--turbo') == (['compass'], '0.95'))
    args = ['compass', '--turbo', '0.95', '--num-processes', '2']
    assert(pop_argument(args, '--num-processes') == (['compass', '--turbo', '0.95'], '2'))
    args = ['compass', '--num-processes', '2', '--turbo', '0.95']
    assert(pop_argument(args, '--num-processes') == (['compass', '--turbo', '0.95'], '2'))
    args = ['compass', '--num-processes', '2']
    assert(pop_argument(args, '--turbo-min-pct') == (['compass', '--num-processes', '2'], ''))


@pytest.mark.slow
def test_compass_resource_manager_get_cell_names():
    r"""
    Tests that the CompassResourceManager correctly finds the cell names.
    """
    compass_args = COMPASS_ARGS
    compass_resource_manager = CompassResourceManager(compass_args=compass_args)
    cell_names = compass_resource_manager.get_cell_names()
    assert(len(cell_names) == 40)  # There should be 40 cells
    assert(cell_names[0] == 'WT-d_S141_L007_R1_001')  # This should be the first cell.


@pytest.mark.slow
def test_compass_resource_manager_get_reaction_ids():
    r"""
    Tests that the CompassResourceManager correctly finds the reaction ids.
    """
    compass_args = COMPASS_ARGS
    compass_resource_manager = CompassResourceManager(compass_args=compass_args)
    reaction_ids = compass_resource_manager.get_reaction_ids()
    assert(len(reaction_ids) == 10211)  # There should be 10211 reactions
    assert(reaction_ids[0] == '10FTHF5GLUtl_pos')  # This should be the first reaction.


@pytest.mark.slow
def test_compass_oracle_shape():
    r"""
    Tests that the CompassOracle returns the correct shape of the underlying matrix.
    Takes 30s on my laptop.
    """
    compass_args = COMPASS_ARGS
    compass_oracle = CompassOracle(compass_args=compass_args)
    matrix_shape = compass_oracle.shape()
    assert(matrix_shape == (40, 10211))  # There should be 40 cells and 10211 reactions.


@pytest.mark.slow
def test_compass_oracle_populate_selected_reactions_file():
    r"""
    Tests that the CompassOracle correctly populates the file passed to
    --selected-reactions-for-each-cell (necessary for running Compass to
    obtain only a subset of the reaction scores.)
    Takes 30s on my laptop.
    """
    compass_args = COMPASS_ARGS
    with tempfile.NamedTemporaryFile("w+") as selected_reactions_file:
        compass_oracle = CompassOracle(compass_args=compass_args)
        compass_oracle._populate_selected_reactions_file(
            [0, 2, 2],
            [0, 0, 1],
            selected_reactions_file.name
        )
        # Check that the written cells and reactions are correct
        print(f"Checking contents of '{selected_reactions_file.name}'")
        with open(selected_reactions_file.name, "r") as file:
            file_content = file.read()
            assert(
                file_content
                == "WT-d_S141_L007_R1_001,10FTHF5GLUtl_pos\n"
                   "WT-1012-c_S135_L007_R1_001,10FTHF5GLUtl_pos,10FTHF5GLUtm_pos\n")


@pytest.mark.slow
def test_compass_oracle_observe_entries():
    r"""
    Tests that the CompassOracle correctly runs Compass under the hood to compute the
    requested reaction scores. This is the key method of the CompassOracle object.
    """
    compass_args = COMPASS_ARGS
    compass_oracle = CompassOracle(compass_args=compass_args)
    vals = compass_oracle.observe_entries(
        [[0, 0], [0, 2], [1, 2]]
    )
    np.testing.assert_almost_equal(vals, np.array([3284.4609895, 3716.8335506, 3620.9471307]))
    vals = compass_oracle.observe_entries(
        [[1, 0]]
    )
    np.testing.assert_almost_equal(vals, np.array([3226.725758506863]))


def test_compass_oracle_cache_observations():
    r"""
    Tests that the caching method of the CompassOracle works as expected.
    """
    with tempfile.NamedTemporaryFile() as output_file:
        CompassOracle._cache_observations(
            [], [], [], output_file.name
        )  # This tests that adding no observations works too.
        CompassOracle._cache_observations(
            [1], [2], [3.14], output_file.name
        )
        CompassOracle._cache_observations(
            [3, 5], [4, 6], [2.71, -1.0], output_file.name
        )
        with open(output_file.name, "r") as f:
            cache_contents = f.read()
            assert(cache_contents == "1 2 3.14\n3 4 2.71\n5 6 -1.0\n")


@pytest.mark.slow
@pytest.mark.parametrize("num_processes", ["1", "2"])
def test_turbo_compass_smoke_test(num_processes):
    r"""
    Just checks that turbo Compass runs end-to-end without crashing in a very
    minimalistic setting.
    Takes 80s on my laptop.
    We use 1 process to turn multi-processing off, which can otherwise
    hide errors. We also try with 2 processes to excercise multiprocessing.
    You can also manually run this from the command line as follows:
    $ compass --media media1 --data compass/Resources/Test\ Data/rsem_tpmTable_2cells.txt \
        --output-dir _test_turbo_compass_smoke_test --num-processes 1 --test-mode \
        --turbo 0.7 --turbo-increments 0.5
    """
    with tempfile.TemporaryDirectory() as output_dir_name:
        print(f"created temporary output directory '{output_dir_name}'")
        compass_args = [
            "compass",
            "--media", "media1",
            "--data", "compass/Resources/Test Data/rsem_tpmTable_2cells.txt",
            "--output-dir", output_dir_name,
            "--num-processes", num_processes,
            "--test-mode",
            "--turbo", "0.7",
            "--turbo-increments", "0.5"
        ]

        with patch.object(sys, 'argv', compass_args):
            compass.main.entry()
        # Try to read the reactions file
        reaction_scores_path = os.path.join(output_dir_name, "reactions.tsv")
        reaction_scores_df = pd.read_csv(reaction_scores_path, sep="\t", index_col=0)
        assert(reaction_scores_df.shape == (106, 2))


@pytest.mark.slow
def test_turbo_compass_smoke_test_2():
    r"""
    Same test as above, but now runs on ALL reactions (does not use --test-mode)
    To avoid this test being slow due to CPLEX, CPLEX is stubbed.
    Takes 75s on my laptop.
    We use 1 process to turn multi-processing off, which can otherwise
    hide errors. Also, we _must_ turn multiprocessing off or else the
    subprocesses won't have CPLEX stubbed and so the test would take forever.
    You can also manually run this from the command line as follows, but do note
    that CPLEX _WON'T BE STUBBED_ so it will take a while:
    $ compass --media media1 --data compass/Resources/Test\ Data/rsem_tpmTable_2cells.txt \
        --output-dir _test_turbo_compass_smoke_test_2 --num-processes 1 \
        --turbo 0.7 --turbo-increments 0.5
    """
    with tempfile.TemporaryDirectory() as output_dir_name:
        print(f"created temporary output directory '{output_dir_name}'")
        compass_args = [
            "compass",
            "--media", "media1",
            "--data", "compass/Resources/Test Data/rsem_tpmTable_2cells.txt",
            "--output-dir", output_dir_name,
            "--num-processes", "1",
            "--turbo", "0.7",
            "--turbo-increments", "0.5"
        ]

        def cplex_solve_stub(self):
            class DummySolution():
                def get_objective_value(self):
                    return 42.0
            self.solution = DummySolution()
        with patch.object(sys, 'argv', compass_args),\
                patch.object(compass.compass.algorithm.cplex.Cplex, 'solve', new=cplex_solve_stub):
            compass.main.entry()
        # Try to read the reactions file
        reaction_scores_path = os.path.join(output_dir_name, "reactions.tsv")
        reaction_scores_df = pd.read_csv(reaction_scores_path, sep="\t", index_col=0)
        assert(reaction_scores_df.shape == (10211, 2))


@pytest.mark.slow
def test_turbo_compass_caching():
    r"""
    Tests that Turbo Compass caching works.
    Takes 3m30s on my laptop.
    """
    with tempfile.TemporaryDirectory() as output_dir_name:
        print(f"created temporary output directory '{output_dir_name}'")
        compass_args = [
            "compass",
            "--media", "media1",
            "--data", "compass/Resources/Test Data/rsem_tpmTable_4cells.txt",
            "--output-dir", output_dir_name,
            "--num-processes", "4",
            "--turbo", "0.7",
            "--turbo-increments", "0.25",
            "--turbo-max-iters", "2"
        ]
        cache_path = os.path.join(output_dir_name, "_tmp", "turbo_compass_cache", "turbo_compass_cache.txt")
        np.random.seed(1)

        def cplex_solve_stub(self):
            class DummySolution():
                def get_objective_value(self):
                    return np.random.uniform(1000.0, 10000.0)
            self.solution = DummySolution()
        with patch.object(sys, 'argv', compass_args),\
                patch.object(compass.compass.algorithm.cplex.Cplex, 'solve', new=cplex_solve_stub):
            # First run
            compass.main.entry()
            # Check that exactly one half of the matrix was queried
            cache = pd.read_csv(cache_path, sep=" ", header=None)
            assert(cache.shape == (2 * 10211, 3))

            # Second run
            compass.main.entry()
            # Check that 3/4 of the matrix was queried
            cache = pd.read_csv(cache_path, sep=" ", header=None)
            assert(cache.shape == (3 * 10211, 3))

            # Third run
            compass.main.entry()
            # Check that all of the matrix was queried
            cache = pd.read_csv(cache_path, sep=" ", header=None)
            assert(cache.shape == (4 * 10211, 3))

            # Fourth run: Even when we observe all the data, Turbo Compass should
            # not crash (i.e. Turbo Compass should be idempotent)
            compass.main.entry()
            # Check that all of the matrix was queried
            cache = pd.read_csv(cache_path, sep=" ", header=None)
            assert(cache.shape == (4 * 10211, 3))


@pytest.mark.slow
def test_turbo_compass_on_40_cells_stubbing_cplex_to_a_rank_1_matrix():
    r"""
    Runs turbo Compass on a 40 cells x 10211 reactions matrix. CPLEX is stubbed
    to produce a rank-1 reaction score matrix.
    Takes 11m on my computer (using 12 processes).
    NOTE: Don't try to run this on the command line since it relies on CPLEX patching!
    (it will take forever if not)
    """
    with tempfile.TemporaryDirectory() as output_dir_name:
        print(f"created temporary output directory '{output_dir_name}'")
        compass_args = [
            "compass",
            "--media", "media1",
            "--data", "compass/Resources/Test Data/rsem_tpmTable_full.txt",
            "--output-dir", output_dir_name,
            "--turbo", "0.95",
            "--turbo-increments", "0.05",
            "--turbo-min-pct", "0.99"
        ]

        def hash_str(s: str, mod: int = 10000) -> int:
            r"""
            Hashes the string modulo 'mod', returning an integer between
            0 and mod-1.
            """
            res = int(hashlib.sha1(s.encode()).hexdigest(), 16) % mod
            return res

        def reaction_score(cell_name: str, reaction_id: str) -> float:
            r"""
            Each cell and reaction has a latent number between 0 and 100.
            Each entry of the reaction score matrix is the product of these
            two latent numbers (i.e. the reaction score matrix has rank 1).
            """
            return np.sqrt(hash_str(cell_name)) * np.sqrt(hash_str(reaction_id))

        def cplex_solve_stub(self):
            class HardcodedSolution():
                def __init__(self, objective_value: float):
                    self.objective_value = objective_value

                def get_objective_value(self) -> float:
                    return self.objective_value
            current_cell_name = global_state.get_current_cell_name()
            current_reaction_id = global_state.get_current_reaction_id()
            self.solution = HardcodedSolution(
                objective_value=reaction_score(
                    current_cell_name,
                    current_reaction_id
                ))
        with patch.object(sys, 'argv', compass_args),\
                patch.object(compass.compass.algorithm.cplex.Cplex, 'solve', new=cplex_solve_stub):
            compass.main.entry()
            # Check that the solution makes sense:
            # (1) Retrieve the imputed solution
            imputed_reaction_scores_df = CompassResourceManager(compass_args).get_reaction_scores()
            # (2) Compute the ground truth solution
            true_reaction_scores_df = imputed_reaction_scores_df.copy()
            for reaction_id in true_reaction_scores_df.index:
                for cell_name in true_reaction_scores_df.columns:
                    true_reaction_scores_df.loc[reaction_id, cell_name] =\
                        reaction_score(cell_name, reaction_id)
            # (3) Compute the spearman correlation for each reaction
            X_true = true_reaction_scores_df.to_numpy().T
            X_imputed = imputed_reaction_scores_df.to_numpy().T
            from turbo_mc.models.cv_matrix_completion_model import compute_spearman_r2s
            cv_spearman_r2s = compute_spearman_r2s(X_true, X_imputed)
            # (4) Check that the median reconstruction quality of the reactions is good.
            # NOTE: Because of caching in Compass' maximize_reaction, constant reactions
            # never hit cplex.solve() and are hardcoded by Compass as 0 (rather than using my
            # low rank approximation). Thus, true_reaction_scores_df is actually wrong for
            # these reactions. They are approximately 30%. To avoid these 30% wrong reactions
            # (that have a SR2 of 0.0), we look at the MEDIAN SR2 in this test, as follows:
            assert(np.sort(cv_spearman_r2s)[len(cv_spearman_r2s) // 2] > 0.95)
