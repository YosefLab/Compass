import tempfile

from compass.global_state import _SelectedReactionsForEachCell


def test_SelectedReactionsForEachCell():
    r"""
    Test that information of which reactions should be computed for each cell is loaded correctly.
    """
    state = _SelectedReactionsForEachCell("tests/data/selected_reactions_for_each_cell.txt")
    assert(state.is_selected("WT-d_S141_L007_R1_001", "10FTHF5GLUtl_pos"))
    assert(not state.is_selected("WT-d_S141_L007_R1_001", "10FTHF5GLUtm_pos"))
    assert(state.is_selected("Ob-DHA-e_S154_L007_R1_001", "10FTHF5GLUtl_pos"))
    assert(state.is_selected("Ob-DHA-e_S154_L007_R1_001", "10FTHF5GLUtm_pos"))
    assert(not state.is_selected("Ob-DHA-b_S132_L007_R1_001", "10FTHF5GLUtl_pos"))
    assert(not state.is_selected("Ob-DHA-b_S132_L007_R1_001", "10FTHF5GLUtm_pos"))
    assert(not state.is_selected("unseen_cell", "10FTHF5GLUtl_pos"))
    assert(not state.is_selected("WT-d_S141_L007_R1_001", "unseen_reaction"))
    assert(not state.is_selected("unseen_cell", "unseen_reaction"))


def test_by_default_all_cells_and_reactions_are_selected():
    r"""
    If no file with selected reactions is provided, all cells and
    reactions should be valid.
    """
    state = _SelectedReactionsForEachCell(None)
    assert(state.is_selected("any_cell", "any_reaction"))


def test_with_empty_file_nothing_is_selected():
    r"""
    If an empty file is supplied, nothing should be computed.
    (This is just a border case I ran into. It basically guarantees that turbo Compass
    is idempotent.)
    """
    with tempfile.NamedTemporaryFile("w+") as file:
        state = _SelectedReactionsForEachCell(file.name)
        assert(not state.is_selected("some_cell", "some_reaction"))
