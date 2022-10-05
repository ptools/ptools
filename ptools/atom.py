"""ptools.atom - Defines classes and function that handle atom and
atom groups."""

# Allows to use type hinting of class within itself
# e.g. BaseAtom.copy(self) -> BaseAtom
from __future__ import annotations

# Pythonn core libraries.
from collections import UserList
import copy
import math

# Scientific libraries.
import numpy as np

# Type-hinting specific import
from typing import Any, Sequence
from numpy.typing import ArrayLike

# PTools imports.
from . import linalg
from . import tables
from .spatial import TransformableObject, TranslatableObject
from .io.formatters.pdb import PDBFormatter

# The Protein Data Bank format for atom coordinates
PDB_FMT = (
    "{record:<6s}{index:5s} "
    "{name:4s}{altloc}{resname:<4s}{chain:s}{resid:>4s}{insertion}   "
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
        resid (int): residue number (read from topology file)
        charge (float): charge
        coords (numpy.ndarray): cartesian coordinates
        meta (dict[str, ()]): metadata dictionnary
    """

    # Atom equality evaluated on these attributes
    _comparison_attributes = ["index", "name", "resname", "resid", "chain", "coords"]

    def __init__(
        self,
        index: int = 0,
        name: str = "XXX",
        resname: str = "XXX",
        chain: str = "X",
        resid: int = 0,
        charge: float = 0.0,
        coords: ArrayLike = np.zeros(3),
        meta: dict = None,
    ):
        super().__init__(coords)
        self._name = name
        self.resname = resname
        self.chain = chain
        self.index = index
        self.resid = resid
        self.charge = charge
        self.meta = meta

    @property
    def element(self):
        """Returns an atom element name (read-only)."""
        return guess_atom_element(self.name)

    @property
    def name(self) -> str:
        """Gets/Sets atom's name."""
        return self._name

    @name.setter
    def name(self, name: str):
        """Name setter simultaneously updates element name."""
        self._name = name

    def __eq__(self, other: object) -> bool:
        """Compares two BaseAtom instances."""
        if not isinstance(other, BaseAtom):
            raise TypeError(f"cannot compare BaseAtom with object of type {type(other)}")
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

    def __init__(self, atom: BaseAtom, serial: int, collection: AtomCollection):
        self.serial = serial
        self.collection = collection
        attrs = (
            "name",
            "resname",
            "chain",
            "index",
            "resid",
            "charge",
            "meta",
            "coords",
        )
        kwargs = {k: getattr(atom, k) for k in attrs}
        super().__init__(**kwargs)

    @property
    def coords(self) -> np.ndarray:
        """Gets atom cartesian coordinates."""
        return self.collection.coords[self.serial].copy()

    @coords.setter
    def coords(self, pos: np.ndarray):
        self.collection.coords[self.serial] = np.array(pos)

    @property
    def mass(self) -> float:
        """Gets/Sets an atom mass."""
        return self.collection.masses[self.serial]

    @mass.setter
    def mass(self, mass: float):
        """Gets/Sets an atom mass."""
        self.collection.masses[self.serial] = mass

    # Override BaseAtom name setter to manually handle masses
    # (name getter override is mandatory).
    @property
    def name(self) -> str:
        """Gets/Sets atom's name."""
        return self._name

    @name.setter
    def name(self, name: str):
        self._name = name
        self.collection.masses[self.serial] = guess_atom_mass(self.element)


class AtomCollection(TransformableObject, UserList):
    """Group of atoms.

    For better performances, atom coordinates are stored into a numpy array.

    Args:
        atoms (list[BaseAtom]): list of atoms
    """

    def __init__(self, atoms: Sequence[BaseAtom] = None):
        if atoms is None or len(atoms) == 0:
            atoms = []
            coords = np.zeros((0, 3))
        else:
            atoms = [Atom(atom, serial, self) for serial, atom in enumerate(atoms)]
            coords = np.array([atom._coords for atom in atoms])

        TransformableObject.__init__(self, coords)
        UserList.__init__(self, atoms)
        self.masses = np.zeros(len(atoms))
        self.guess_masses()

    def __repr__(self) -> str:
        """String representation."""
        modulename = self.__module__
        classname = self.__class__.__name__
        return f"<{modulename}.{classname} with {len(self)} atoms>"

    def __add__(self, other: AtomCollection) -> AtomCollection:
        """Concatenates two RigidBody instances."""
        output = super().__add__(other.copy())
        output.coords = np.concatenate((self.coords, other.coords), axis=0)
        return output

    def __iadd__(self, other: AtomCollection) -> AtomCollection:
        return self.__add__(other)

    def guess_masses(self):
        """Guesses atom masses and store them."""
        self.masses = np.array([guess_atom_mass(atom.element) for atom in self])

    def copy(self) -> AtomCollection:
        """Returns a copy of the current collection."""
        return self.__class__(self)

    def size(self) -> int:
        """Gets the number of atoms in the collection.

        Alias for len(AtomCollection).
        """
        return len(self)

    def center_to_origin(
        self, origin: ArrayLike = np.zeros(3), use_weights: bool = False
    ):
        """Centers AtomCollection on `origin`."""
        if not use_weights:
            super().center_to_origin(origin)
        else:
            self.translate(np.array(origin) - self.center_of_mass())

    def center_of_mass(self) -> np.ndarray:
        """Returns the center of mass (barycenter)."""
        return linalg.center_of_mass(self.coords, self.masses)

    def inertia_tensor(self, weights=None):
        """Returns the inertia tensors of a set of atoms."""
        if weights is None:
            weights = self.masses
        return linalg.inertia_tensor(self.coords, weights)

    def principal_axes(self, sort: bool = True) -> np.ndarray:
        """Returns an AtomCollection principal axes.

        Args:
            sort (bool): sort axes by importance
        """
        return linalg.principal_axes(self.inertia_tensor(), sort)

    def radius_of_gyration(self) -> float:
        """Returns the isometric radius of gyration (atom mass is not taken
        into account)."""
        centered = self.coords - self.centroid()
        rgyr2 = np.sum(centered**2) / len(self)
        return math.sqrt(rgyr2)

    def topdb(self) -> str:
        """Returns a string representing the AtomCollection in PDB format."""
        return "\n".join(atom.topdb() for atom in self)

    def writepdb(self, path: str):
        """Writes the AtomCollection to a PDB formatted file."""
        with open(path, "wt", encoding="utf-8") as f:
            print(self.topdb(), file=f)

    def set_chain(self, chain: str):
        """Sets all atom chain property."""
        for atom in self:
            atom.chain = chain

    def select_atom_type(self, atom_type: str) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired atom type."""
        return AtomCollection(atoms=[atom for atom in self if atom.name == atom_type])

    def select_atom_types(self, atom_types: list[str]) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired atom types."""
        return AtomCollection(
            atoms=[atom for atom in self if atom.name.strip() in atom_types]
        )

    def select_residue_range(self, start: int, end: int) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired which residue is within the range."""
        return AtomCollection(
            atoms=[atom for atom in self if start <= atom.resid <= end]
        )

    def select_chain(self, chain_id: str) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired chain."""
        return AtomCollection(atoms=[atom for atom in self if atom.chain == chain_id])


def guess_atom_element(atom_name: str) -> str:
    """Returns the atom element based on its name.

    Basically returns the first non-numeric character in atom_name.
    """
    for char in atom_name:
        if char.isalpha():
            return char
    return "X"


def guess_atom_mass(element: str) -> float:
    """Returns the atom mass based on the element name."""
    return tables.masses.get(element, 1.0)
