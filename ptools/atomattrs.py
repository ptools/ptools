"""ptools.atom - Defines classes and function that handle atom and
atom groups."""

from __future__ import annotations

# Python core libraries.
import copy

# Type-hinting specific import
from typing import Any, Optional

# Scientific libraries.
import numpy as np

# Other third-party libraries.
from attrs import define, field

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
    coordinates: array3d = field(
        factory=lambda: array3d((0, 0, 0)), converter=lambda x: array3d(x)
    )
    meta: dict[str, Any] = field(factory=dict)
    element: str = ""
    mass: float = -1.0  # -1.0 means "not set"
    radius: float = -1.0  # -1.0 means "not set"

    def __attrs_post_init__(self):
        """Post-initialization method."""
        if not self.element:
            self.element = guess_atom_element(self.name)
        if self.mass < 0:
            self.mass = guess_atom_mass(self.element)
        if self.radius < 0:
            self.radius = guess_atom_radius(self.element)

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
    return tables.atomic_masses.get(element, 1.0)


def guess_atom_element(name: str) -> str:
    """Returns the atom element based on its name.

    Basically returns the first non-numeric character in atom_name.
    """
    for char in name:
        if char.isalpha():
            return char
    return "X"


def guess_atom_radius(element: str) -> float:
    """Returns the atom radius based on the element name.

    Returns 0.0 if the element is not found in the table.
    """
    return tables.atomic_radii.get(element, 0.0)
