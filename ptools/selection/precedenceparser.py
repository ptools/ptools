"""Implements an operator-precedence parser for the selection language."""

from abc import ABC, abstractmethod
from typing import Protocol, Union

from .errors import (
    UnexpectedTokenError,
    UnexpectedTrailingTokenError,
    UnknownTokenError
)


class Operator(Protocol):

    token: str
    precedence: int
    operands: int

    def eval(self, *args):
        ...


class LogicOperator(Operator):
    """Base class for an logic operator (e.g. 'and', 'or', 'not', etc.)."""


class LeafOperator(Operator):
    """Base class for a leaf operator."""

    precedence = 1
    operands = 1


class EvaluatorBase(ABC):
    """Base class for the selection language parsing."""

    def __init__(self, logic_operators: list[LogicOperator]):
        self.binary_operators = {op.token: op for op in logic_operators if op.operands == 2}
        self.unary_operators = {op.token: op for op in logic_operators if op.operands == 1}
        self.leaf_operators: dict[str, Operator] = {}
        self.tokens: list[str] = []
        self.cursor: int = 0

    def register_leaf_operator(self, operator: LeafOperator):
        self.leaf_operators[operator.token] = operator

    def is_keyword(self, token: str) -> bool:
        return (
            token in self.leaf_operators or
            token in self.binary_operators or
            token in self.unary_operators
        )

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

    def _expect(self, token: str):
        if self._next() == token:
            self._consume()
        else:
            raise UnexpectedTokenError(token, self._next())  # type: ignore

    def _expect_end(self):
        try:
            self._expect(None)
        except UnexpectedTokenError as excp:
            raise UnexpectedTrailingTokenError(self._next()) from excp

    def _eval_node(self, operator, *args):
        return operator.eval(*args)

    def _eval_leaf(self, token):
        values = []
        while self._has_next() and not self.is_keyword(self._next()):
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
        # mypy is actually quite dumb on this one: not able to properly infer that arguments
        # are cannot be None when they reach each particular instruction.
        t = self._parse_subexpression()
        while (
            self._next() in self.binary_operators and
                self.binary_operators[self._next()].precedence >= precedence  # type: ignore
        ):
            op = self.binary_operators.get(self._consume())
            t1 = self.parse_expression(op.precedence + 1)  # type: ignore
            t = self._eval_node(op, t, t1)
        return t

    def _parse_subexpression(self):
        if self._next() in self.unary_operators:
            op = self.unary_operators[self._consume()]
            t = self.parse_expression(op.precedence)
            return self._eval_node(op, t)

        # At this point, the token should be a leaf.
        if self._next() and self._next() in self.leaf_operators:
            t = self._eval_leaf(self._consume())
            return t

        # At this point, the keyword has not been understood.
        raise UnknownTokenError(self._next())
        return None
