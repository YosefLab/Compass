r"""
This module contains global state. It is used for e.g. determining whether
to skip the computation of a reaction for a cell when --selected-reactions-for-each-cell
is provided.
"""
import os
from typing import Optional

_current_cell_name = -1
_current_reaction_id = ""
# This will hold an object of the class SelectedReactionsForEachCell that
# will allow us to know
_selected_reactions_for_each_cell = None


class _SelectedReactionsForEachCell():
    def __init__(self, path: Optional[str]):
        r"""
        :param path: Path to file containing the reactions for each cell.
            The file contains one line per cell of interest. Each line is
            comma-separated. It starts with the name of the cell, and is
            followed by the reaction ids of interest.
            E.g. 'cell42,BUP2_pos,BUP2_neg,DHPM2_pos' is a valid line.
            If no path is provided, then no (cells, reaction) pairs will
            be excluded by this method.
        """
        if path is None:
            reactions_for_cells = None
        else:
            reactions_for_cells = {}
            if not os.path.exists(path):
                raise Exception(f"cannot find selected reactions for each cell: file {path}")

            with open(path) as f:
                for line in f:
                    # Must be careful to remove newlines from end of file.
                    tokens = [token.rstrip("\n") for token in line.split(',')]
                    cell_name, reactions_for_this_cell = tokens[0], tokens[1:]
                    reactions_for_cells[cell_name] = reactions_for_this_cell
        self.reactions_for_cells = reactions_for_cells

    def is_selected(self, cell_name: str, reaction_id: str) -> bool:
        r"""
        Whether the reaction_id is selected for cell_name.
        """
        if self.reactions_for_cells is None:
            # No file was supplied: everything is selected by default.
            return True
        return reaction_id in self.reactions_for_cells.get(cell_name, [])


def init_selected_reactions_for_each_cell(path: Optional[str]) -> None:
    global _selected_reactions_for_each_cell
    _selected_reactions_for_each_cell = _SelectedReactionsForEachCell(path)


def set_current_cell_name(cell_name: str) -> None:
    global _current_cell_name
    _current_cell_name = cell_name


def set_current_reaction_id(reaction_id: str) -> None:
    global _current_reaction_id
    _current_reaction_id = reaction_id


def get_current_cell_name() -> str:
    return _current_cell_name


def get_current_reaction_id() -> str:
    return _current_reaction_id


def current_reaction_is_selected_for_current_cell() -> bool:
    return _selected_reactions_for_each_cell.is_selected(_current_cell_name, _current_reaction_id)
