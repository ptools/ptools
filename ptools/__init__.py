"""Top-level package for PTools"""

__version__ = "0.2.0"

from .io import read_pdb as read_pdb
from .io import write_pdb as write_pdb

from .attract import AttractRigidBody
from .particlecollection import ParticleCollection
from .rigidbody import RigidBody


from . import (
    atomattrs,
    attract,
    forcefield,
    heligeom,
    io,
    pairlist,
    particlecollection,
    reduce,
    rigidbody,
    selection,
    superpose,
    tables,
)


__all__ = [
    "atomattrs",
    "attract",
    "forcefield",
    "heligeom",
    "io",
    "pairlist",
    "particlecollection",
    "reduce",
    "rigidbody",
    "selection",
    "superpose",
    "tables",
    "AttractRigidBody",
    "RigidBody",
]
