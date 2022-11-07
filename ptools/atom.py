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
from .spatial import TranslatableObject
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
class BaseAtom(TranslatableObject):
    """Base class for an Atom.

    Args:
        index (int): index (read from topology file)
        name (str): atom name
        resname (str): residue name
        chain (str): chain identifier
        residue_index (int): residue number (read from topology file)
        charge (float): charge
        coords (numpy.ndarray): cartesian coordinates
        meta (dict[str, ()]): metadata dictionnary
    """

    # Atom equality evaluated on these attributes
    _comparison_attributes = ["index", "name", "residue_name", "residue_index", "chain", "coords"]

    def __init__(
        self,
        index: int = 0,
        name: str = "XXX",
        residue_name: str = "XXX",
        chain: str = "X",
        residue_index: int = 0,
        charge: float = 0.0,
        meta: dict = None,
        coords: ArrayLike = np.zeros(3),
    ):
        super().__init__(coords)
        self.name = name
        self.index = index
        self.residue_index = residue_index
        self.residue_name = residue_name
        self.chain = chain
        self.charge = charge
        self.meta = meta

    @property
    def element(self):
        """Returns an atom element name (read-only)."""
        return self.guess_element(self.name)

    def __eq__(self, other: object) -> bool:
        """Compares two BaseAtom instances."""
        if not isinstance(other, BaseAtom):
            raise TypeError(
                f"cannot compare BaseAtom with object of type {type(other)}"
            )
        for key in self._comparison_attributes:
            if key != "coords":
                if getattr(self, key) != getattr(other, key):
                    return False
            elif np.abs(self.coords - other.coords).sum() > 1e-6:
                return False
        return True

    def __repr__(self) -> str:
        """BaseAtom string representation."""
        attrs = {}
        for key, value in sorted(self.__dict__.items()):
            if key.startswith("_") and not key.startswith("__"):
                key = key[1:]
            attrs[key] = value
        modulename = self.__module__
        classname = self.__class__.__name__
        desc = ", ".join(f"{key}={value!r}" for key, value in attrs.items())
        return f"<{modulename}.{classname}({desc})>"

    def copy(self) -> BaseAtom:
        """Returns a copy of the current atom."""
        obj = copy.copy(self)
        obj.coords = obj.coords.copy()
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
            "coords",
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
