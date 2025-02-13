"""Top-level package for PTools"""

from loguru import logger

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
from .attract import AttractDockingParameters, AttractRigidBody
from .io import read_pdb as read_pdb
from .io import write_pdb as write_pdb
from .particlecollection import ParticleCollection
from .rigidbody import RigidBody
from .transform import move as move

__version__ = "0.2.1"


__all__ = [
    "atomattrs",
    "attract",
    "forcefield",
    "heligeom",
    "io",
    "move",
    "pairlist",
    "particlecollection",
    "reduce",
    "rigidbody",
    "selection",
    "superpose",
    "tables",
    "read_pdb",
    "write_pdb",
    "AttractDockingParameters",
    "AttractRigidBody",
    "RigidBody",
    "ParticleCollection",
]

logger.disable("ptools")


class RigidBodyFactory:
    @staticmethod
    def from_pdb(path: str) -> RigidBody:
        return RigidBody(read_pdb(path))
