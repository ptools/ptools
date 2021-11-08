"""Defines classes and functions used for the multicopy feature."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from ptools.forcefield import AttractForceField1, ForceFieldBase
from ptools.rigidbody import RigidBodyBase

from .atom import AtomCollection
from .io import pdb


class McopRigidPDBError(Exception):
    """Raised when an error occurs while reading a McopRigid PDB."""


@dataclass
class Mcop(RigidBodyBase):
    """Container for multiple copies of the same monomer.

    Inherits from RigidBodyBase which ensures that it implements the `from_pdb` method.
    """

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

    def copy(self) -> Mcop:
        """Returns a copy of itself."""
        output = Mcop()
        output.copies = [atoms.copy() for atoms in self.copies]
        return output

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
class McopRigid(RigidBodyBase):
    """Mcop RigidBody."""

    core: AtomCollection = None
    regions: Sequence[Mcop] = field(default_factory=list)
    weights: Sequence[np.ndarray] = field(default_factory=list)

    def copy(self) -> McopRigid:
        """Returns a copy of itself."""
        output = McopRigid()
        output.core = self.core.copy()
        output.regions = [region.copy() for region in self.regions]
        output.weights = [weights.copy() for weights in self.weights]

    def add_region(self, region: Mcop):
        """Add a region to the internal list of regions."""
        self.regions.append(region)
        self.weights.append(np.zeros(len(region)))

    def center(self, origin: np.ndarray = np.zeros(3)):
        """Center core and regions on `origin`."""
        self.core.center(origin)
        for region in self.regions:
            region.center(origin)

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
        self.core = models[keys[0][0]]

        # Then should come regions.
        if len(keys) < 2:
            raise McopRigidPDBError("no region found")
        for i, model in enumerate(models.values()):
            if i > 0:
                self.regions.append(model)


class LigandNotInitializedError(ValueError):
    """Raised where trying to set receptor but ligand has not been set previously."""


@dataclass
class McopForceField(ForceFieldBase):

    cutoff: float = field(init=False, default=10.0)
    _ff_core: AttractForceField1 = field(init=False, default=None)
    _centered_ligand: McopRigid = field(init=False, default=None)
    _moved_ligand: McopRigid = field(init=False, default=None)
    _receptor: McopRigid = field(init=False, default=None)

    def energy(self) -> float:
        return 0.0

    def update(self):
        self.energy()

    def set_receptor(self, receptor: McopRigid):
        if self._moved_ligand is None:
            raise LigandNotInitializedError("Trying to set receptor while ligand has not been set")
        self._receptor = receptor.copy()
        self._receptor.center(origin=self._moved_ligand.core.center_of_mass())
        self._init_energies()

    def set_ligand(self, ligand: McopRigid):
        self._moved_ligand = ligand.copy()
        self._centered_ligand = ligand.copy()
        self._centered_ligand.center(origin=ligand.core.center_of_mass())


    def _init_energies(self):
        self._ff_core = AttractForceField1(self._receptor.core, self._moved_ligand.core, self.cutoff)

    def non_bonded_energy(self) -> float:
        energy_region = 0.0
        energy_core = 0.0

        print(self._ff_core.non_bonded_energy())


