"""ptools.atom - Defines classes and function that handle atom and
atom groups."""

from __future__ import annotations

# Python core libraries.
import copy

# Scientific libraries.
import numpy as np

# Other third-party libraries.
from attrs import define, field

# Type-hinting specific import
from typing import Any, Optional

# PTools imports.
from . import tables
from .array3d import array3d


# pylint: disable=R0902,R0913
# A lot of instant attributes... Is it really an issue?
@define(slots=False, eq=False)
class AtomAttrs:
    """Stores atom properties."""

    name: str = "XXX"
    index: int = 0
    residue_name: str = "XXX"
    residue_index: int = 0
    chain: str = "X"
    charge: float = 0.0
    coordinates: array3d = field(
        factory=lambda: array3d((0, 0, 0)), converter=lambda x: array3d(x)
    )
    meta: dict[str, Any] = field(factory=dict)
    element: Optional[str] = None
    mass: Optional[float] = None

    def __attrs_post_init__(self):
        """Post-initialization method."""
        if self.element is None:
            self.element = guess_atom_element(self.name)
        if self.mass is None:
            self.mass = guess_atom_mass(self.element)

    def __eq__(self, other: object) -> bool:
        """Tests for equality.

        ``coordinates`` equality is check using ``numpy.allclose``.
        """
        if not isinstance(other, AtomAttrs):
            raise NotImplementedError

        attrs = self.__dict__.keys() - ("coordinates",)
        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False

        return np.allclose(self.coordinates, other.coordinates)

    def copy(self) -> AtomAttrs:
        """Returns a copy of the current atom."""
        obj = copy.deepcopy(self)
        return obj

    def guess_mass(self):
        """Returns the atom mass based on the element name."""
        self.mass = guess_atom_mass(self.element)

    def guess_element(self):
        """Returns the atom element based on its name.

        Basically returns the first non-numeric character in atom_name.
        """
        self.element = guess_atom_element(self.name)


def guess_atom_mass(element: str) -> float:
    """Returns the atom mass based on the element name."""
    return tables.masses.get(element, 1.0)


def guess_atom_element(name: str) -> str:
    """Returns the atom element based on its name.

    Basically returns the first non-numeric character in atom_name.
    """
    for char in name:
        if char.isalpha():
            return char
    return "X"
