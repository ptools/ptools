"""Defines classes and functions used for the multicopy feature."""

import numpy as np

from dataclasses import dataclass, field
from typing import Sequence

from .atom import Atom, AtomCollection
from .io import read_pdb


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
        models = read_pdb(path)
        assert isinstance(models, list)
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
