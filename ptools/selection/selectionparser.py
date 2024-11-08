"""ptools.selection - Selection language and parsing."""
from __future__ import annotations

import re
from typing import Any, Callable, Protocol, TYPE_CHECKING, Union

import numpy as np

from .precedenceparser import (
    PrecedenceClimbingEvaluator,
    LeafOperator,
    LogicOperator,
)

if TYPE_CHECKING:
    from .particlecollection import ParticleCollection

Numeric = Union[int, float]

# ===========================================================================
#
# Logic Operators
#
# ===========================================================================


class AndOperator:
    """Selection 'and' operator."""

    token = "and"
    precedence = 3
    operands = 2

    def eval(self, *args):
        assert len(args) == self.operands
        left, right = args
        indices = np.intersect1d(left.serial, right.serial)
        return left[indices]


class OrOperator:
    """Selection 'or' operator."""

    token = "or"
    precedence = 3
    operands = 2

    def eval(self, *args):
        assert len(args) == self.operands
        return args[0] + args[1]


class NotOperator:
    """Selection 'not' operator. """

    token = "not"
    precedence = 5
    operands = 1

    def eval(self, *args):
        assert len(args) == self.operands
        notsel = args[0]
        indices = np.setdiff1d(notsel.parent.serial, notsel.serial)
        return notsel.parent[indices]


# ===========================================================================
#
# Property Selection operators.
#
# For selecting atoms based on their properties.
# For example, 'resid 1 2 3' selects atoms with residue index 1, 2, or 3.
#
# ===========================================================================


class PropertySelectionOperator:
    """Base class for selections based on atom properties."""

    precedence: int
    operands: int
    token: str
    attr: str

    def __init__(self, token: str, attr: str):
        ...

    def eval(self, atoms: "ParticleCollection", values: list):
        ...


class IntegerPropertySelection:
    """Base class for selections based on integer attribute (e.g. atom or residue index)."""

    precedence = 1
    operands = 1

    def __init__(self, token: str, attr: str):
        self.token = token
        self.attr = attr

    def eval(self, atoms: "ParticleCollection", values: list[Numeric]):
        ops = {
            ">": np.greater,
            "<": np.less,
            ">=": np.greater_equal,
            "<=": np.less_equal,
            "==": np.equal,
            "!=": np.not_equal,
        }

        if values[0] in ops:
            return self.select_from_operator(atoms, ops[values[0]], int(values[1]))

        # Range selection in fashion attr <number> to <number>.
        if len(values) == 3 and values[1] == "to":
            start, end = int(values[0]), int(values[2])
            return self.select_from_range(atoms, start, end)

        values = [int(resid) for resid in values]
        return self.select_from_values(atoms, values)

    def select_from_operator(self, atoms: "ParticleCollection", op: Callable, value: Numeric):
        """Selection from a single operator."""
        indices = np.where(op(atoms.atom_properties.get(self.attr).values, value))[0]
        return atoms[indices]

    def select_from_values(self, atoms: "ParticleCollection", value: list[Numeric]):
        """Selection from a single or multiple values."""
        indices = np.where(np.isin(atoms.atom_properties.get(self.attr).values, value))[0]
        return atoms[indices]

    def select_from_range(self, atoms: "ParticleCollection", start: Numeric, end: Numeric):
        """Selection from a range of values."""
        indices = np.where(
            np.logical_and(
                atoms.atom_properties.get(self.attr).values >= start,
                atoms.atom_properties.get(self.attr).values <= end,
            )
        )[0]
        return atoms[indices]


class StringPropertySelection:
    """Base class for selections based on string attribute (e.g. atom or residue name)."""

    operands = 1
    precedence = 1

    def __init__(self, token: str, attr: str):
        self.token = token
        self.attr = attr

    def eval(self, atoms: "ParticleCollection", values: list[str]):
        indices = np.where(np.isin(atoms.atom_properties.get(self.attr).values, values))[0]
        return atoms[indices]


class BoolPropertySelection:
    """Base class for selections based on boolean attribute (e.g. hetero)."""

    operands = 1
    precedence = 1

    def __init__(self, token: str, attr: str):
        self.token = token
        self.attr = attr

    def eval(self, atoms: "ParticleCollection", values: list[str]):
        indices = np.where(atoms.atom_properties.get(self.attr).values == True)[0]
        return atoms[indices]


# ===========================================================================
#
# Keyword Selection
#
# ===========================================================================
class WaterSelection:
    """Selection operator for water molecules"""

    operands = 0
    precedence = 1
    token = "water"

    def eval(self, atoms: "ParticleCollection", values: list[Any]):
        assert len(values) == self.operands
        water_residues = ["HOH", "WAT", "TIP3", "TIP4", "TIP5"]
        indices = np.where(np.isin(atoms.atom_properties.get("residue_names").values, water_residues))[0]
        return atoms[indices]


# ===========================================================================
#
# Actual Ptools Selection Parser
#
# ===========================================================================
class SelectionParser(PrecedenceClimbingEvaluator):
    def __init__(self, atoms: "ParticleCollection"):
        super().__init__(logic_operators=[AndOperator(), OrOperator(), NotOperator()])
        self.atoms = atoms

        # Pattern for tokenizing the selection string, i.e. separating
        # arithmetic operators from other elements (e.g. 'resid<5' -> 'resid < 5').
        pattern = r'\w+|<=|>=|==|!=|[+\-*/=<>()]'
        self._tokenize_regex = re.compile(pattern)

        # Dynamically creates selection operators based on the ParticleCollection
        # atom properties.
        # Importantly, the selection string uses the singular form of particle
        # properties, e.g. 'names' -> 'name CA'.
        for prop in self.atoms.atom_properties:
            if np.issubdtype(prop.values.dtype, np.number):
                self.register_leaf_operator(IntegerPropertySelection(prop.singular, prop.plural))
            elif np.issubdtype(prop.values.dtype, np.bool_):
                self.register_leaf_operator(BoolPropertySelection(prop.singular, prop.plural))
            else:
                self.register_leaf_operator(StringPropertySelection(prop.singular, prop.plural))

        # Aliases for 'residue_index' and 'residue_name'.
        self.leaf_operators["resid"] = self.leaf_operators["residue_index"]
        self.leaf_operators["resname"] = self.leaf_operators["residue_name"]

        # Registers keyword operators.
        self.register_leaf_operator(WaterSelection())

    def parse(self, selection_str: str):
        """Parses and evaluates the selection string."""
        self.tokens = self._tokenize_regex.findall(selection_str)
        return self.evaluate()

    def _eval_leaf(self, token):
        values = super()._eval_leaf(token)
        operator = self.leaf_operators[token]
        return operator.eval(self.atoms, values)


def select(selection_str: str, atoms: ParticleCollection):
    """Selection function."""
    parser = SelectionParser(atoms)
    selection_str = selection_str.replace(":", " to ")
    return parser.parse(selection_str)
