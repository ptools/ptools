"""ptools.selection - Selection language and parsing."""

from typing import Any

import numpy as np


from .precedenceparser import (
    PrecedenceClimbingEvaluator,
    LeafOperator,
    LogicOperator,
)


# ===========================================================================
#
# Logic Operators
#
# ===========================================================================


class AndOperator(LogicOperator):
    """Selection 'and' operator."""

    token = "and"
    precedence = 3
    operands = 2

    def eval(self, *args):
        assert len(args) == self.operands
        left, right = args
        indices = np.intersect1d(left.serial, right.serial)
        return left[indices]


class OrOperator(LogicOperator):
    """Selection 'or' operator."""

    token = "or"
    precedence = 3
    operands = 2

    def eval(self, *args):
        assert len(args) == self.operands
        return args[0] + args[1]


class NotOperator(LogicOperator):
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
# Selection operators
#
# ===========================================================================

class SelectionOperator(LeafOperator):

    def __init__(self, token: str, attr: str):
        """Initializes the selection operator.

        Attrs:
            token: the token used in the selection string.
            attr: the attribute name in the ParticleCollection.atom_properties dictionary.
              Typically, `token` is singular and `attr` is plural.
        """
        self.token = token
        self.attr = attr


class IntegerAttributeSelection(SelectionOperator):
    """Base class for selections based on integer attribute (e.g. atom or residue index)."""

    def eval(self, *args):
        assert len(args) == 3
        attr, atoms, values = args

        # Range selection in fashion attr <number>:<number>.
        if len(values) == 1 and ":" in values[0]:
            start, end = values[0].split(":")
            return self.select_from_range(atoms, attr, int(start), int(end))

        # Range selection in fashion attr <number> to <number>.
        if len(values) == 3 and values[1] == "to":
            start, end = int(values[0]), int(values[2])
            return self.select_from_range(atoms, attr, start, end)

        values = [int(resid) for resid in values]
        return self.select_from_values(atoms, attr, values)

    @classmethod
    def select_from_values(cls, atoms: "ParticleCollection", attr: str, value: list[Any]):
        """Selection from a single or multiple values."""
        indices = np.where(np.isin(atoms.atom_properties.get(attr).values, value))[0]
        return atoms[indices]

    @classmethod
    def select_from_range(cls, atoms, attr, start, end):
        """Selection from a range of values."""
        indices = np.where(
            np.logical_and(
                atoms.atom_properties.get(attr).values >= start,
                atoms.atom_properties.get(attr).values <= end,
            )
        )[0]
        return atoms[indices]


class StringAttributeSelection(SelectionOperator):
    """Base class for selections based on integer attribute (e.g. atom or residue name)."""

    # def eval(self, token, atoms, values):
    def eval(self, *args):
        assert len(args) == 3
        attr, atoms, values = args
        indices = np.where(np.isin(atoms.atom_properties.get(attr).values, values))[0]
        return atoms[indices]


class SelectionParser(PrecedenceClimbingEvaluator):
    def __init__(self, atoms: "ParticleCollection"):
        super().__init__(logic_operators=[AndOperator(), OrOperator(), NotOperator()])
        self.atoms = atoms

        # Dynamically creates selection operators based on the ParticleCollection
        # atom properties.
        # Importantly, the selection string uses the singular form of particle
        # properties, e.g. 'names' -> 'name CA'.

        for prop in self.atoms.atom_properties:
            if np.issubdtype(prop.values.dtype, np.number):
                self.leaf_operators[prop.singular] = IntegerAttributeSelection(prop.singular, prop.plural)
            else:
                self.leaf_operators[prop.singular] = StringAttributeSelection(prop.singular, prop.plural)

        # Aliases for 'residue_index' and 'residue_name'.
        self.leaf_operators["resid"] = self.leaf_operators["residue_index"]
        self.leaf_operators["resname"] = self.leaf_operators["residue_name"]

    def parse(self, selection_str: str):
        self.tokens = selection_str.split()
        return self.evaluate()

    def _eval_leaf(self, token):
        values = super()._eval_leaf(token)
        operator = self.leaf_operators[token]
        return operator.eval(operator.attr, self.atoms, values)


def select(selection_str: str, atoms: "ParticleCollection"):
    """Selection function."""
    parser = SelectionParser(atoms)
    return parser.parse(selection_str)
