"""ptools.atom - Defines classes and function that handle atom and
atom groups."""

# Allows to use type hinting of class within itself
# e.g. BaseAtom.copy(self) -> BaseAtom
from __future__ import annotations

# Python core libraries.
import copy
from dataclasses import dataclass, field

# Scientific libraries.
import numpy as np

# Type-hinting specific import
from typing import Any, TYPE_CHECKING

from ptools.array3d import array3d
from ._typing import ArrayLike

# PTools imports.
from . import tables
from .io.formatters.pdb import PDBFormatter



if TYPE_CHECKING:
    from atomcollection import AtomCollection



# The Protein Data Bank format for atom coordinates
PDB_FMT = (
    "{record:<6s}{index:5s} "
    "{name:4s}{altloc}{resname:<4s}{chain:s}{residue_index:>4s}{insertion}   "
    "{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          "
    "{element:>2s}"
)


# pylint: disable=R0902,R0913
# A lot of instant attributes... Is it really an issue?
@dataclass
class BaseAtom:
    """Base class for an Atom."""

    name: str = "XXX"
    index: int = 0
    residue_name: str = "XXX"
    residue_index: int = 0
    chain: str = "X"
    charge: float = 0.0
    coordinates: array3d = field(default_factory=lambda: array3d((0, 0, 0)))
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        self.coordinates = array3d(self.coordinates)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BaseAtom):
            raise NotImplementedError

        for attr, value in self.__dict__.items():
            if attr != "coordinates":
                if value != getattr(other, attr):
                    return False
            elif not np.allclose(self.coordinates, other.coordinates):
                return False

        return True

    @property
    def coords(self) -> array3d:
        return self.coordinates

    @coords.setter
    def coords(self, value: ArrayLike):
        self.coordinates = array3d(value)

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



@dataclass(init=False)
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