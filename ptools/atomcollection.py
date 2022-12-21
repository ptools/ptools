from __future__ import annotations

from collections import UserList
import itertools
from typing import Any, Callable, Iterable, Iterator, Optional, TypeVar

import numpy as np

from .array3d import array3d
from .atomattrs import AtomAttrs
from ._typing import ArrayLike


class Atom(AtomAttrs):
    """Atom that belongs to a group of atoms.

    Its coordinates are a weak reference to the AtomCollection coordinate
    array. Atom knows its positions in the AtomCollection and can therefore
    find its coordinates into the AtomCollection coordinate array.

    An Atom initialization can only be done from a AtomAttrs.

    Args:
        serial (int): atom index in the collection (can be different from
            `Atom.index` attribute).
        atom (AtomAttrs): atom properties stored in AtomAttrs.
        collection (AtomCollection): reference to the collection the atom
            belongs to.

    """

    def __init__(self, atom: AtomAttrs, serial: int, collection: "AtomCollection"):
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
        return self.collection.coordinates[self.serial].copy()

    @coords.setter
    def coords(self, pos: ArrayLike):
        self.collection.coordinates[self.serial] = array3d(pos)

    @property
    def mass(self) -> float:
        """Gets/Sets an atom mass."""
        return self.collection.masses[self.serial]

    @mass.setter
    def mass(self, mass: float):
        """Gets/Sets an atom mass."""
        self.collection.masses[self.serial] = mass

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, self.__class__):
            err = f"cannot compare {self.__class__.__qualname__} with object of type {type(other)}"
            raise TypeError(err)

        if not np.allclose(self.coordinates, other.coordinates):
            return False

        attrs = self.__dict__.keys() - ("coordinates", "collection")
        return all(getattr(self, attr) == getattr(other, attr) for attr in attrs)


AtomCollectionType = TypeVar("AtomCollectionType", bound="AtomCollection")

class AtomCollection(UserList):
    """Group of atoms.

    For better performances, atom coordinates are stored into a numpy array.

    Args:
        atoms (list[AtomAttrs]): list of atoms
    """

    def __init__(self, atoms: Optional[Iterable[AtomAttrs]] = None):
        if atoms is None:
            atoms = []
            coords = array3d(np.zeros((0, 3)))

        atoms = [Atom(atom, serial, self) for serial, atom in enumerate(atoms)]
        if atoms:
            coords = array3d([atom.coordinates for atom in atoms])

        self.coordinates = coords
        UserList.__init__(self, atoms)
        self.masses = np.zeros(len(atoms))
        self.guess_masses()

    def __eq__(self, other: object) -> bool:
        """Checks for equality."""
        if not isinstance(other, self.__class__):
            err = f"cannot compare {self.__class__.__qualname__} with object of type {type(other)}"
            raise TypeError(err)

        if not np.allclose(self.coordinates, other.coordinates):
            return False
        return all(lhs == rhs for (lhs, rhs) in zip(self, other))

    def __repr__(self) -> str:
        """String representation."""
        modulename = self.__module__
        classname = self.__class__.__name__
        return f"<{modulename}.{classname} with {len(self)} atoms>"

    def __add__(self: AtomCollectionType, other: Iterable[AtomAttrs]) -> AtomCollectionType:
        """Concatenates two RigidBody instances."""
        if not isinstance(other, AtomCollection):
            other = AtomCollection(other)
        output = super().__add__(other.copy())
        lhs = np.asarray(self.coordinates)
        rhs = np.asarray(other.coordinates)
        output.coordinates = np.concatenate((lhs, rhs), axis=0)
        return output

    def __iadd__(self: AtomCollectionType, other: Iterable[AtomAttrs]) -> AtomCollectionType:
        return self.__add__(other)

    def guess_masses(self):
        """Guesses atom masses and store them."""
        self.masses = np.array([Atom.guess_mass(atom.element) for atom in self])

    def copy(self: AtomCollectionType) -> AtomCollectionType:
        """Returns a copy of the current collection."""
        return self.__class__(self)

    def size(self) -> int:
        """Gets the number of atoms in the collection.

        Alias for len(AtomCollection).
        """
        return len(self)

    def set_chain(self, chain: str):
        """Sets all atom chain property."""
        for atom in self:
            atom.chain = chain

    def groupby(self: AtomCollectionType, key: Callable) -> dict[Any, AtomCollectionType]:
        data = sorted(self, key=key)
        grouped = itertools.groupby(data, key=key)
        return {key: self.__class__(list(group)) for key, group in grouped}

    def select_atom_type(self: AtomCollectionType, atom_type: str) -> AtomCollectionType:
        """Returns a sub-collection made of atoms with desired atom type."""
        return self.__class__(atoms=[atom for atom in self if atom.name == atom_type])

    def select_atom_types(self: AtomCollectionType, atom_types: list[str]) -> AtomCollectionType:
        """Returns a sub-collection made of atoms with desired atom types."""
        return self.__class__(atoms=[atom for atom in self if atom.name in atom_types])

    def select_residue_range(self: AtomCollectionType, start: int, end: int) -> AtomCollectionType:
        """Returns a sub-collection made of atoms with desired which residue is within the range."""
        return self.__class__(
            atoms=[atom for atom in self if start <= atom.residue_index <= end]
        )

    def select_chain(self: AtomCollectionType, chain_id: str) -> AtomCollectionType:
        """Returns a sub-collection made of atoms with desired chain."""
        return self.__class__(atoms=[atom for atom in self if atom.chain == chain_id])

    def iter_atoms(self) -> Iterator[Atom]:
        """Iterate over the collection's atoms."""
        return iter(self)

    def iter_residues(self: AtomCollectionType) -> Iterator[AtomCollectionType]:
        by_residue = self.groupby(lambda atom: (atom.residue_index, atom.chain))
        return iter(by_residue.values())

    def iter_chains(self: AtomCollectionType) -> Iterator[AtomCollectionType]:
        by_chain = self.groupby(lambda atom: atom.chain)
        return iter(by_chain.values())
