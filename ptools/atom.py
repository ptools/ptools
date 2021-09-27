"""ptools.atom - Defines classes and function that handle atom and
atom groups."""

# Allows to use type hinting of class within itself
# e.g. BaseAtom.copy(self) -> BaseAtom
from __future__ import annotations

import collections.abc
import copy
import math
from typing import Iterator, Sequence, Union

import numpy as np

import ptools.spatial
from ptools.spatial import SpatialObject
from . import tables


# The Protein Data Bank format for atom coordinates
PDB_FMT = (
    "{record:<6s}{index:5s} "
    "{name:4s}{altloc}{resname:<4s}{chain:s}{resid:>4s}{insertion}   "
    "{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          "
    "{element:>2s}"
)


# pylint: disable=R0902,R0913
# A lot of instant attributes... Is it really an issue?
class BaseAtom(SpatialObject):
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

    def __init__(
        self,
        index: int = 0,
        name: str = "XXX",
        resname: str = "XXX",
        chain: str = "X",
        resid: int = 0,
        charge: float = 0.0,
        coords: np.ndarray = np.zeros(3),
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
    
    def __eq__(self, other: BaseAtom) -> bool:
        """Compares two BaseAtom instances."""
        for key, value in self.__dict__.items():
            if key.startswith("_") and not key.startswith("__"):
                key = key[1:]
            if key != "coords" and value != getattr(other, key):
                return False
            if key == "coords":
                if np.abs(self.coords - other.coords).sum() > 1e-6:
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

        rec = "ATOM"
        altlocchar = " "
        insertion = " "
        occ = 1.0
        bfactor = 0.0
        element = self.element

        chain = " " if not self.chain else self.chain
        namebuf = f"{self.name:>4s}" if len(self.name) > 2 else self.name.center(4)
        indexbuf = f"{self.index:5d}" if self.index < 100000 else f"{self.index:05x}"
        residbuf = f"{self.resid:4d}" if self.resid < 10000 else f"{self.resid:04x}"

        x, y, z = self.coords

        return PDB_FMT.format(
            record=rec,
            index=indexbuf,
            name=namebuf,
            altloc=altlocchar,
            resname=self.resname,
            chain=chain[0],
            resid=residbuf,
            insertion=insertion[0],
            x=x,
            y=y,
            z=z,
            occupancy=occ,
            bfactor=bfactor,
            element=element,
        )


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
    
    def __eq__(self, other: BaseAtom) -> bool:
        """Compares two Atom instances."""
        for key, value in self.__dict__.items():
            if key.startswith("_") and not key.startswith("__"):
                key = key[1:]
            if key not in ("coords", "collection") and value != getattr(other, key):
                return False
            if key == "coords":
                if np.abs(self.coords - other.coords).sum() > 1e-6:
                    return False
        return True


class AtomCollection(SpatialObject, collections.abc.Sequence):
    """Group of atoms.

    For better performances, atom coordinates are stored into a numpy array.

    Args:
        atoms (list[BaseAtom]): list of atoms
    """

    def __init__(self, atoms: Sequence[BaseAtom] = None):
        if atoms is None:
            atoms = []
        self.atoms = [Atom(atom, serial, self) for serial, atom in enumerate(atoms)]
        if self.atoms:
            coords = np.array([atom._coords for atom in self.atoms])
        else:
            coords = np.zeros((0, 3))
        super().__init__(coords)
        self.masses = np.zeros(len(self.atoms))
        self.guess_masses()

    def __contains__(self, __x: Atom) -> bool:
        return __x in self.atoms

    def __len__(self) -> int:
        """Gets the number of atoms in the collection."""
        return len(self.atoms)

    def __repr__(self) -> str:
        """String representation."""
        modulename = self.__module__
        classname = self.__class__.__name__
        return f"<{modulename}.{classname} with {len(self)} atoms>"

    def __iter__(self) -> Iterator[Atom]:
        """Iterates over the collection atoms."""
        return iter(self.atoms)

    def __getitem__(self, serial: Union[int, slice]) -> Union[Atom, AtomCollection]:
        """Accesses an atom by its serial number (which the internal index
        starting at 0)."""
        if isinstance(serial, slice):
            return self.__class__(atoms=self.atoms[serial])
        return self.atoms[serial]

    def __add__(self, other: AtomCollection) -> AtomCollection:
        """Concatenates two RigidBody instances."""
        output = self.copy()
        other = other.copy()
        output.atoms += list(other)
        output.coords = np.concatenate((output.coords, other.coords), axis=0)
        return output

    def guess_masses(self):
        """Guesses atom masses and store them."""
        self.masses = np.array([guess_atom_mass(atom.element) for atom in self])

    def copy(self) -> AtomCollection:
        """Returns a copy of the current collection."""
        return self.__class__(self.atoms)

    def size(self) -> int:
        """Gets the number of atoms in the collection.

        Alias for len(AtomCollection).
        """
        return len(self)

    def center(self) -> np.ndarray:
        """Returns the isobarycenter (geometric center) of a collection of
        atoms."""
        return self.centroid()

    def center_of_mass(self) -> np.ndarray:
        """Returns the center of mass (barycenter)."""
        return ptools.spatial.center_of_mass(self.coords, self.masses)

    def tensor_of_inertia(self, weights=None, method: str = "accurate"):
        """Returns the inertia tensors of a set of atoms.

        Args:
            method (str): "fast" or "accurate"
                The "fast" method does not take into account atom masses.
        """
        weights = self.masses if method == "accurate" else None
        return super().tensor_of_inertia(weights, method)

    def principal_axes(self, sort: bool = True, method: str = "accurate") -> np.ndarray:
        """Returns an AtomCollection principal axes.

        Args:
            sort (bool): sort axes by importance
            method (str): "fast" or "accurate"
                The "fast" method does not take into account atom masses when
                calculating the tensor of inertia, resulting in a less
                accurate result.
                The "fast" method is probably sufficient to calculate axes of
                inertia.
        """
        return ptools.spatial.principal_axes(self.tensor_of_inertia(method), sort)

    def radius_of_gyration(self) -> float:
        """Returns the isometric radius of gyration (atom mass is not taken
        into account)."""
        centered = self.coords - self.center()
        rgyr2 = np.sum(centered ** 2) / len(self)
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
        for atom in self.atoms:
            atom.chain = chain

    def select_atom_type(self, atom_type: str) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired atom type."""
        return self.__class__(atoms=[atom for atom in self if atom.name == atom_type])

    def select_residue_range(self, start: int, end: int) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired which residue is within the range."""
        return self.__class__(
            atoms=[atom for atom in self if start <= atom.resid <= end]
        )

    def select_chain(self, chain_id: str) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired chain."""
        return self.__class__(atoms=[atom for atom in self if atom.chain == chain_id])


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
