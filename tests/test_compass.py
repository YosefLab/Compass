import os
import numpy as np
import pandas as pd
import pytest
import sys
import tempfile
from unittest.mock import patch

import compass.main


@pytest.mark.slow
@pytest.mark.parametrize("num_processes", ["1", "2"])
def test_compass_smoke_test(num_processes):
    r"""
    Just checks that Compass runs end-to-end without crashing in a very
    minimalistic setting.
    Takes 54s on my laptop.
    We use 1 process to turn multi-processing off, which can otherwise
    hide errors. We also try with 2 processes to excercise multiprocessing.
    You can also manually run this from the command line as follows:
    $ compass --media media1 --data compass/Resources/Test\ Data/rsem_tpmTable_2cells.txt \
        --output-dir _test_compass_smoke_test --num-processes 1 --test-mode
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
        ]
        with patch.object(sys, 'argv', compass_args):
            compass.main.entry()
        # Try to read the reactions file
        reaction_scores_path = os.path.join(output_dir_name, "reactions.tsv")
        reaction_scores_df = pd.read_csv(reaction_scores_path, sep="\t", index_col=0)
        assert(reaction_scores_df.shape == (106, 2))


@pytest.mark.slow
def test_compass_smoke_test_2():
    r"""
    Same test as above, but now runs on ALL reactions (does not use --test-mode)
    To avoid this test being slow due to CPLEX, CPLEX is stubbed.
    Takes 30s on my laptop.
    We use 1 process to turn multi-processing off, which can otherwise
    hide errors. Also, we _must_ turn multiprocessing off or else the
    subprocesses won't have CPLEX stubbed and so the test would take forever.
    You can also manually run this from the command line as follows, but do note
    that CPLEX _WON'T BE STUBBED_ so it will take a while:
    $ compass --media media1 --data compass/Resources/Test\ Data/rsem_tpmTable_2cells.txt \
        --output-dir _test_compass_smoke_test_2 --num-processes 1
    """
    with tempfile.TemporaryDirectory() as output_dir_name:
        print(f"created temporary output directory '{output_dir_name}'")
        compass_args = [
            "compass",
            "--media", "media1",
            "--data", "compass/Resources/Test Data/rsem_tpmTable_2cells.txt",
            "--output-dir", output_dir_name,
            "--num-processes", "1"
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
def test_compass_select_reactions_argument():
    r"""
    Runs Compass with the '--select_reactions' argument in a minimalistic
    setting, checking that the outputted reaction scores are correct.
    Takes 10s on my laptop.
    """
    with tempfile.TemporaryDirectory() as output_dir_name:
        print(f"created temporary output directory '{output_dir_name}'")
        compass_args = [
            "compass",
            "--media", "media1",
            "--data", "compass/Resources/Test Data/rsem_tpmTable_2cells.txt",
            "--output-dir", output_dir_name,
            "--single-sample", "1",
            "--select_reactions", "tests/data/select_reactions.txt"]
        with patch.object(sys, 'argv', compass_args):
            compass.main.entry()
        # Check that contents of outputted reaction scores are correct
        reaction_scores_path = os.path.join(output_dir_name, "_tmp", "reactions.txt")
        reaction_scores_df = pd.read_csv(reaction_scores_path, sep="\t", index_col=0)
        assert(list(reaction_scores_df.index) == ["10FTHF5GLUtl_pos"])
        assert(list(reaction_scores_df.columns) == ["Ob-DHA-e_S154_L007_R1_001"])
        np.testing.assert_almost_equal(reaction_scores_df.to_numpy(), np.array([[3226.7257585068633]]))


@pytest.mark.slow
def test_compass_select_reactions_for_each_cell_argument():
    r"""
    Runs Compass with the '--selected-reactions-for-each-cell' argument in a minimalistic
    setting, checking that the outputted reaction scores are correct.
    Takes 25s on my laptop.
    """
    with tempfile.TemporaryDirectory() as output_dir_name:
        print(f"created temporary output directory '{output_dir_name}'")
        compass_args = [
            "compass",
            "--media", "media1",
            "--data", "compass/Resources/Test Data/rsem_tpmTable_4cells.txt",
            "--selected-reactions-for-each-cell", "tests/data/selected_reactions_for_each_cell.txt",
            "--output-dir", output_dir_name
        ]
        with patch.object(sys, 'argv', compass_args):
            compass.main.entry()
        # Check that contents of outputted reaction scores are correct
        reaction_scores_path = os.path.join(output_dir_name, "reactions.tsv")
        reaction_scores_df = pd.read_csv(reaction_scores_path, sep="\t", index_col=0)
        assert(len(reaction_scores_df.index) == 10211)
        # I will only check the top 3 x 4 rows of the matrix.
        assert(list(reaction_scores_df.index)[:3] == [
            "10FTHF5GLUtl_pos",
            "10FTHF5GLUtm_pos",
            "10FTHF6GLUtl_pos"
        ])
        assert(list(reaction_scores_df.columns) == [
            "WT-d_S141_L007_R1_001",
            "Ob-DHA-e_S154_L007_R1_001",
            "WT-1012-c_S135_L007_R1_001",
            "Ob-DHA-b_S132_L007_R1_001"
        ])
        np.testing.assert_almost_equal(
            reaction_scores_df.iloc[:3, :],
            np.array(
                [[3284.4609894639148, 3226.7257585068633, np.nan, np.nan],
                    [np.nan, 4142.645986483619, np.nan, np.nan],
                    [np.nan, np.nan, np.nan, np.nan]]))
