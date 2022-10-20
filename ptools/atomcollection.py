from __future__ import annotations

from collections import UserList
import itertools
import math
from typing import Any, Callable, Iterator, Sequence

import numpy as np

from .atom import Atom, BaseAtom
from .spatial import TransformableObject
from . import linalg

from ._typing import ArrayLike, FilePath


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

    def writepdb(self, path: FilePath):
        """Writes the AtomCollection to a PDB formatted file."""
        with open(path, "wt", encoding="utf-8") as f:
            print(self.topdb(), file=f)

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
        return self.__class__(
            atoms=[atom for atom in self if atom.name in atom_types]
        )

    def select_residue_range(self, start: int, end: int) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired which residue is within the range."""
        return self.__class__(
            atoms=[atom for atom in self if start <= atom.resid <= end]
        )

    def select_chain(self, chain_id: str) -> AtomCollection:
        """Returns a sub-collection made of atoms with desired chain."""
        return self.__class__(atoms=[atom for atom in self if atom.chain == chain_id])

    def iter_atoms(self) -> Iterator[Atom]:
        """Iterate over the collection's atoms."""
        return iter(self)

    def iter_residues(self) -> Iterator[AtomCollection]:
        by_residue = self.groupby(lambda atom: (atom.resid, atom.chain))
        return iter(by_residue.values())

    def iter_chains(self) -> Iterator[AtomCollection]:
        by_chain = self.groupby(lambda atom: atom.chain)
        return iter(by_chain.values())
