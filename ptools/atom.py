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


@dataclass
class ChainProperties:
    name: str = "X"

    def copy(self) -> ChainProperties:
        return copy.deepcopy(self)


@dataclass
class ResidueProperties:
    """Stores residue properties."""
    name: str = "XXX"
    index: int = 0

    def copy(self) -> ResidueProperties:
        return copy.deepcopy(self)

@dataclass
class AtomProperties:
    """Stores atom properties."""
    name: str = "XXX"
    index: int = 0
    residue: ResidueProperties = field(default_factory = ResidueProperties)
    chain: ChainProperties = field(default_factory = ChainProperties)
    charge: float = 0.0
    meta: dict[str, Any] = field(default_factory=dict)

    def copy(self) -> AtomProperties:
        return copy.deepcopy(self)


# pylint: disable=R0902,R0913
# A lot of instant attributes... Is it really an issue?
class BaseAtom(TranslatableObject):
    """Base class for an Atom."""

    def __init__(
        self,
        index: int = 0,
        name: str = "XXX",
        residue_name: str = "XXX",
        chain: str = "X",
        residue_index: int = 0,
        charge: float = 0.0,
        meta: dict[str, Any] = None,
        coords: ArrayLike = np.zeros(3),
    ):
        super().__init__(coords)
        self.properties = AtomProperties()
        self.properties.name = name
        self.properties.index = index
        self.properties.charge = charge
        self.properties.meta = meta
        self.properties.residue.index = residue_index
        self.properties.residue.name = residue_name
        self.properties.chain.name = chain


    def __getattribute__(self, __name: str) -> Any:
        if __name == "properties":
            return super().__getattribute__(__name)

        if hasattr(self, "properties"):
            attrs = super().__getattribute__("properties")

            if __name == "residue_index":
                return attrs.residue.index
            if __name == "residue_name":
                return attrs.residue.name
            if __name == "chain":
                return attrs.chain.name

            if __name in attrs.__dict__:
                return getattr(attrs, __name)

        return super().__getattribute__(__name)


    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name == "properties":
            return super().__setattr__(__name, __value)

        if hasattr(self, "properties"):
            attrs = super().__getattribute__("properties")

            if __name == "residue_index":
                return setattr(attrs.residue, "index", __value)
            if __name == "residue_name":
                return setattr(attrs.residue, "name", __value)
            if __name == "chain":
                return setattr(attrs.chain, "name", __value)

            if __name in attrs.__dict__:
                return setattr(attrs, __name, __value)


        return super().__setattr__(__name, __value)


    @property
    def element(self):
        """Returns an atom element name (read-only)."""
        return self.guess_element(self.name)

    def __eq__(self, other: object) -> bool:
        """Compares two BaseAtom instances."""
        if not isinstance(other, BaseAtom):
            err = f"cannot compare BaseAtom with object of type {type(other)}"
            raise TypeError(err)

        if self.properties != other.properties:
            return False

        return np.allclose(self.coords, other.coords)

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
        obj = copy.deepcopy(self)
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
