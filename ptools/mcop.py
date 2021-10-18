"""Defines classes and functions used for the multicopy feature."""

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from .atom import AtomCollection
from .io import pdb


class McopRigidPDBError(Exception):
    """Raised when an error occurs while reading a McopRigid PDB."""


@dataclass
class Mcop:
    """Container for multiple copies of the same monomer."""

    copies: Sequence[AtomCollection] = field(default_factory=list)

    def add_copy(self, kopy: AtomCollection):
        """Append a copy to the list of copies."""
        self.copies.append(kopy)

    def size(self):
        """Returns the number of copies.

        Alias:
            len(mcop)
        """
        return len(self.copies)

    def __len__(self):
        """Returns the number of copies.

        Alias:
            mcop.size()
        """
        return self.size()

    def __getitem__(self, key: int) -> AtomCollection:
        """Returns the copy at key `key`."""
        return self.copies[key]

    def clear(self):
        """Clears the internal list of copies."""
        self.copies.clear()

    def read_pdb(self, path: str):
        """Initializes models from data read in a PDB file."""
        models = pdb.read_pdb(path)
        for mod in models:
            self.add_copy(mod)

    def attract_euler_rotate(self, phi: float, ssi: float, rot: float):
        """Rotates copies with Attract Euler protocol."""
        for kopy in self.copies:
            kopy.attract_euler_rotate(phi, ssi, rot)


@dataclass
class McopRigid:
    """Mcop RigidBody."""

    core: AtomCollection = None
    regions: Sequence[Mcop] = field(default_factory=list)
    weights: Sequence[np.ndarray] = field(default_factory=list)

    def add_region(self, region: Mcop):
        """Add a region to the internal list of regions."""
        self.regions.append(region)
        self.weights.append(np.zeros(len(region)))

    def attract_euler_rotate(self, phi: float, ssi: float, rot: float):
        """Rotates core and copies with Attract Euler protocol."""
        self.core.attract_euler_rotate(phi, ssi, rot)
        for region in self.regions:
            region.attract_euler_rotate(phi, ssi, rot)

    def read_pdb(self, path: str):
        """Initializez McopRigid from data read in PDB file."""
        models = pdb.read_pdb(path, as_dict=True)
        keys = [key.split() for key in models.keys()]

        # Core region should come first.
        if len(keys[0]) > 1:
            raise McopRigidPDBError(f"expecting core region first (found f{keys[0]})")

        # Then should come regions.
        if len(keys) < 2:
            raise McopRigidPDBError("no region found")
