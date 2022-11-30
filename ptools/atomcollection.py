from __future__ import annotations

from collections import UserList
import itertools
import math
from typing import Any, Callable, Iterable, Iterator

import numpy as np


from .array3d import array3d
from .atomattrs import AtomAttrs
from .spatial import SupportsTransformation
from . import linalg
from . import measure
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
        return super().__eq__(other)


class AtomCollection(SupportsTransformation, UserList):
    """Group of atoms.

    For better performances, atom coordinates are stored into a numpy array.

    Args:
        atoms (list[AtomAttrs]): list of atoms
    """

    def __init__(self, atoms: Iterable[AtomAttrs] = None):
        if atoms is None:
            atoms = []
            coords = array3d(np.zeros((0, 3)))

        atoms = [Atom(atom, serial, self) for serial, atom in enumerate(atoms)]
        if atoms:
            coords = array3d([atom.coordinates for atom in atoms])

        SupportsTransformation.__init__(self, coords)
        UserList.__init__(self, atoms)
        self.masses = np.zeros(len(atoms))
        self.guess_masses()

    def __repr__(self) -> str:
        """String representation."""
        modulename = self.__module__
        classname = self.__class__.__name__
        return f"<{modulename}.{classname} with {len(self)} atoms>"

    def __add__(self, other: Iterable[AtomAttrs]) -> AtomCollection:
        """Concatenates two RigidBody instances."""
        if not isinstance(other, AtomCollection):
            other = AtomCollection(other)
        output = super().__add__(other.copy())
        output.coordinates = np.concatenate((self.coordinates, other.coordinates), axis=0)
        return output

    def __iadd__(self, other: Iterable[AtomAttrs]) -> AtomCollection:
        return self.__add__(other)

    def guess_masses(self):
        """Guesses atom masses and store them."""
        self.masses = np.array([Atom.guess_mass(atom.element) for atom in self])

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
            self.translate(np.array(origin) - measure.center_of_mass(self))

    def set_chain(self, chain: str):
        """Sets all atom chain property."""
        for atom in self:
            atom.chain = chain

    def groupby(self, key: Callable) -> dict[Any, AtomCollection]:
        data = sorted(self, key=key)
        grouped = itertools.groupby(data, key=key)
        return {key: self.__class__(list(group)) for key, group in grouped}

    def select_atom_type(self, atom_type: str) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired atom type."""
        return self.__class__(atoms=[atom for atom in self if atom.name == atom_type])

    def select_atom_types(self, atom_types: list[str]) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired atom types."""
        return self.__class__(atoms=[atom for atom in self if atom.name in atom_types])

    def select_residue_range(self, start: int, end: int) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired which residue is within the range."""
        return self.__class__(
            atoms=[atom for atom in self if start <= atom.residue_index <= end]
        )

    def select_chain(self, chain_id: str) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired chain."""
        return self.__class__(atoms=[atom for atom in self if atom.chain == chain_id])

    def iter_atoms(self) -> Iterator[Atom]:
        """Iterate over the collection's atoms."""
        return iter(self)

    def iter_residues(self) -> Iterator[AtomCollection]:
        by_residue = self.groupby(lambda atom: (atom.residue_index, atom.chain))
        return iter(by_residue.values())

    def iter_chains(self) -> Iterator[AtomCollection]:
        by_chain = self.groupby(lambda atom: atom.chain)
        return iter(by_chain.values())
