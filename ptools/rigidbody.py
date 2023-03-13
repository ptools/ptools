"""ptools.rigidbody - Defines the RigidBody class and children."""

from __future__ import annotations

# Type hinting specific imports
from collections.abc import Sequence
from typing import Optional, Type, TypeVar
from ._typing import FilePath

# Scientific libraries.
import numpy as np

# PTools.
from .atomattrs import AtomAttrs
from .io.pdb import read_pdb as io_read_pdb
from .particlecollection import ParticleCollection


RigidBodyType = TypeVar("RigidBodyType", bound="RigidBody")

class RigidBody(ParticleCollection):
    """RigidBody is an ParticleCollection that can be initialized from a PDB file.

    It can be initialized from a file.
    """

    def __init__(self, atoms: Optional[Sequence[AtomAttrs]] = None):
        if isinstance(atoms, str):
            raise TypeError(
                "RigidBody class can not longer be instantiated from a path. "
                "Use RigidBody.from_pdb instead."
            )
        ParticleCollection.__init__(self, atoms)

    @classmethod
    def from_pdb(cls: Type[RigidBodyType], path: FilePath) -> RigidBodyType:
        atom_container = io_read_pdb(path)
        assert isinstance(atom_container, ParticleCollection)
        return cls.from_properties(atom_container.atom_properties)
