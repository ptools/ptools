"""ptools.atom - Defines classes and function that handle atom and
atom groups."""

# Allows to use type hinting of class within itself
# e.g. BaseAtom.copy(self) -> BaseAtom
from __future__ import annotations

# Python core libraries.
import copy
from attrs import define, field

# Scientific libraries.
import numpy as np

# Type-hinting specific import
from typing import Any, TYPE_CHECKING
from ._typing import ArrayLike

# PTools imports.
from . import tables
from .array3d import array3d, validate_array3d_dimensions
from .io.formatters.pdb import PDBFormatter


if TYPE_CHECKING:
    from atomcollection import AtomCollection


# pylint: disable=R0902,R0913
# A lot of instant attributes... Is it really an issue?
@define(slots=False, eq=False)
class BaseAtom:
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

    def __eq__(self, other: object) -> bool:
        """Tests for equality.

        ``coordinates`` equality is check using ``numpy.allclose``.
        """
        if not isinstance(other, BaseAtom):
            raise NotImplementedError

        attrs = self.__dict__.keys() - ("coordinates",)
        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False

        return np.allclose(self.coordinates, other.coordinates)

    @property
    def element(self):
        """Returns an atom element name (read-only)."""
        return self.guess_element(self.name)

    def copy(self) -> BaseAtom:
        """Returns a copy of the current atom."""
        obj = copy.deepcopy(self)
        return obj

    def topdb(self) -> str:
        """Returns the atom's description in PDB format."""
        return PDBFormatter.format_atom(self)

    @classmethod
    def guess_mass(cls, element: str) -> float:
        """Returns the atom mass based on the element name."""
        return tables.masses.get(element, 1.0)

    @classmethod
    def guess_element(cls, name: str) -> str:
        """Returns the atom element based on its name.

        Basically returns the first non-numeric character in atom_name.
        """
        for char in name:
            if char.isalpha():
                return char
        return "X"


class Atom(BaseAtom):
    """Atom that belongs to a group of atoms.

    Its coordinates are a weak reference to the AtomCollection coordinate
    array. Atom knows its positions in the AtomCollection and can therefore
    find its coordinates into the AtomCollection coordinate array.

    An Atom initialization can only be done from a BaseAtom.

    Args:
        serial (int): atom index in the collection (can be different from
            `Atom.index` attribute).
        atom (BaseAtom): atom properties stored in BaseAtom.
        collection (AtomCollection): reference to the collection the atom
            belongs to.

    """

    def __init__(self, atom: BaseAtom, serial: int, collection: "AtomCollection"):
        self.serial = serial
        self.collection = collection
        attrs = (
            "name",
            "residue_name",
            "chain",
            "index",
            "residue_index",
            "charge",
            "meta",
            "coordinates",
        )
        kwargs = {k: getattr(atom, k) for k in attrs}
        super().__init__(**kwargs)

    @property
    def coords(self) -> array3d:
        """Gets atom cartesian coordinates."""
        return self.collection.coords[self.serial].copy()

    @coords.setter
    def coords(self, pos: ArrayLike):
        self.collection.coords[self.serial] = array3d(pos)

    @property
    def mass(self) -> float:
        """Gets/Sets an atom mass."""
        return self.collection.masses[self.serial]

    @mass.setter
    def mass(self, mass: float):
        """Gets/Sets an atom mass."""
        self.collection.masses[self.serial] = mass

    def __eq__(self, other: object) -> bool:
        return super().__eq__(other)
