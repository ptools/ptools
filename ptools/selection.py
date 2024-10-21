"""ptools.selection - Selection language and parsing."""


from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Union

import numpy as np

from .particlecollection import Particle, ParticleCollection


# Translation from selection language tokens to ParticleCollection attributes.
TOKEN_TO_ATTR = {
    "index": "indices",
    "name": "names",
    "resid": "residue_indices",
    "resname": "residue_names",
    "chain": "chains",
}


class UnexpectedTokenError(Exception):
    """Error raised when a token was expected and found another one."""

    def __init__(self, expected: str, found: str):
        self.expected = expected
        self.found = found
        message = f"Expected '{self.expected}', got '{self.found}'"
        super().__init__(message)


class UnexpectedTrailingTokenError(Exception):
    def __init__(self, token: str):
        self.token = token
        message = f"Unexpected token at the end of the selection string: '{self.token}'"
        super().__init__(message)


class MisformattedExpressionError(Exception):
    message_fmt = "Misformatted expression at: '{}'"

    def __init__(self, token: str):
        self.token = token
        super().__init__(self.message_fmt.format(self.token))


class UnknownTokenError(MisformattedExpressionError):
    message_fmt = "Unknown token: '{}'"


class MisformattedRangeError(MisformattedExpressionError):
    pass


class OperatorBase(ABC):
    """Base class for an operator."""

    token: str
    precedence: int
    operands: int = 2

    def is_unary(self) -> bool:
        return self.operands == 1

    def is_binary(self) -> bool:
        return self.operands == 2

    @abstractmethod
    def eval(self, *args):
        """Evaluation: should be override by children classes."""


class AndOperator(OperatorBase):
    """Selection 'and' operator."""
    token = "and"
    precedence = 3

    def eval(self, *args: tuple(ParticleCollection, ParticleCollection)):
        left, right = args
        indices = np.intersect1d(left.serial, right.serial)
        return left[indices]


class OrOperator(OperatorBase):
    """Selection 'or' operator."""
    token = "or"
    precedence = 3

    def eval(self, *args: tuple(ParticleCollection, ParticleCollection)):
        left, right = args
        return left + right


class NotOperator(OperatorBase):
    """Selection 'not' operator."""
    token = "not"
    precedence = 5
    operands = 1

    def eval(self, notsel: ParticleCollection):
        indices = np.setdiff1d(notsel.parent.serial, notsel.serial)
        return notsel.parent[indices]


OPERATORS = [AndOperator(), OrOperator(), NotOperator()]
BINARY_OPERATORS = {op.token: op for op in OPERATORS if op.is_binary()}
UNARY_OPERATORS = {op.token: op for op in OPERATORS if op.is_unary()}

SELECTION_TOKENS = {}


class SelectionMeta(type):
    """Metaclass for Selection classes.

    Only purpose is to automatically register itself to SELECTION_TOKENS.
    """

    def __init__(cls, class_name, bases, attrs):
        type.__init__(type, class_name, bases, attrs)
        try:
            SELECTION_TOKENS[attrs["token"]] = cls
        except KeyError:
            # A SelectionBase child class may not have a "token" attribute
            # if it is supposed to be a parent class for actual selection
            # classes itself.
            pass


class SelectionBase(metaclass=SelectionMeta):
    """Base class for a Selection."""

    token: str
    precedence: int

    def eval(self, atoms: ParticleCollection, values: Any):
        """Evaluation: should be override by children classes."""
        raise NotImplementedError("This method should be overriden in subclasses")


class IntegerAttributeSelection(SelectionBase):
    """Base class for selections based on integer attribute (e.g. atom or residue index)."""

    def eval(self, atoms, values):
        if attr := TOKEN_TO_ATTR.get(self.token) is None:
            raise UnknownTokenError(self.token)
        attr = TOKEN_TO_ATTR[self.token]

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
    def select_from_values(cls, atoms: ParticleCollection, attr: str, value: list[Any]):
        indices = np.where(np.isin(atoms.atom_properties.get(attr).values, value))[0]
        return atoms[indices]

    @classmethod
    def select_from_range(cls, atoms, attr, start, end):
        indices = np.where(
            np.logical_and(
                atoms.atom_properties.get(attr).values >= start,
                atoms.atom_properties.get(attr).values <= end,
            )
        )[0]
        return atoms[indices]


class StringAttributeSelection(SelectionBase):
    """Base class for selections based on integer attribute (e.g. atom or residue name)."""

    def eval(self, atoms, values):
        if attr := TOKEN_TO_ATTR.get(self.token) is None:
            raise UnknownTokenError(self.token)
        attr = TOKEN_TO_ATTR[self.token]
        indices = np.where(np.isin(atoms.atom_properties.get(attr).values, values))[0]
        return atoms[indices]


class ResidueIndexSelection(IntegerAttributeSelection):
    """Selection by residue index."""

    token = "resid"
    precedence = 1


class ParticleIndexSelection(IntegerAttributeSelection):
    """Selection by atom index."""

    token = "index"
    precedence = 1


class ChainSelection(StringAttributeSelection):
    """Selection by chain identifier."""

    token = "chain"
    precedence = 1


class ParticleNameSelection(StringAttributeSelection):
    """Selection by atom name."""

    token = "name"
    precedence = 1


class ResidueNameSelection(StringAttributeSelection):
    """Selection by residue name."""

    token = "resname"
    precedence = 1


def binary(token):
    """Returns the operator if `token` is a binary operator, else None."""
    return BINARY_OPERATORS.get(token)


def unary(token):
    """Returns the operator if `token` is an unary operator, else None."""
    return UNARY_OPERATORS.get(token)


def is_keyword(token):
    """Returns True if a token is a reserved keyword."""
    return (
        token in SELECTION_TOKENS
        or token in BINARY_OPERATORS
        or token in UNARY_OPERATORS
    )


class EvaluatorBase(ABC):
    """Base class for the selection language parsing."""

    def __init__(self):
        self.tokens = []
        self.cursor = 0

    def _next(self) -> Union[str, None]:
        if self.cursor < len(self.tokens):
            return self.tokens[self.cursor]
        return None

    def _has_next(self):
        return self._next() is not None

    def _consume(self):
        token = self._next()
        self.cursor += 1
        return token

    def _error(self, msg=None):
        raise MisformattedExpressionError(msg)

    def _expect(self, token: str):
        if self._next() == token:
            self._consume()
        else:
            raise UnexpectedTokenError(token, self._next())

    def _expect_end(self):
        try:
            self._expect(None)
        except UnexpectedTokenError as excp:
            raise UnexpectedTrailingTokenError(self._next()) from excp

    def _eval_node(self, operator, *args):
        return operator.eval(*args)

    def _eval_leaf(self, token):
        values = []
        while self._has_next() and not is_keyword(self._next()):
            values.append(self._consume())
        return values

    @abstractmethod
    def evaluate(self):
        pass


class PrecedenceClimbingEvaluator(EvaluatorBase):
    """Implements the Precedence Climbing algorithm to parse the selection language."""

    def evaluate(self):
        tree = self.parse_expression(0)
        self._expect_end()
        return tree

    def parse_expression(self, precedence: int):
        t = self._parse_subexpression()
        while (
            self._next() in BINARY_OPERATORS
            and binary(self._next()).precedence >= precedence
        ):
            op = binary(self._consume())
            t1 = self.parse_expression(op.precedence + 1)
            t = self._eval_node(op, t, t1)
        return t

    def _parse_subexpression(self):
        if self._next() in UNARY_OPERATORS:
            op = unary(self._consume())
            t = self.parse_expression(op.precedence)
            return self._eval_node(op, t)

        if self._next() and self._next() in SELECTION_TOKENS:
            t = self._eval_leaf(self._consume())
            return t

        # At this point, the keyword has not been understood
        raise UnknownTokenError(self._next())
        return None


class SelectionParser(PrecedenceClimbingEvaluator):
    def __init__(self, atoms: ParticleCollection = None):
        super().__init__()
        self.atoms = atoms

    def parse(self, selection_str: str, atoms: ParticleCollection = None):
        if atoms is not None:
            self.atoms = atoms
        # if self.atoms is None:
        #     raise ValueError("atoms member should be initialized")
        self.tokens = selection_str.split()
        return self.evaluate()

    def _eval_leaf(self, token):
        values = super()._eval_leaf(token)
        return SELECTION_TOKENS[token]().eval(self.atoms, values)


def select(selection_str: str, atoms: ParticleCollection = None):
    """Selection function."""
    parser = SelectionParser(atoms)
    return parser.parse(selection_str)
