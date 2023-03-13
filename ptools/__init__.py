"""Top-level package for PTools"""

__author__ = """Benoist LAURENT"""
__email__ = "benoistlaurent@gmail.com"
__version__ = "0.1.0"

from .io import read_pdb as read_pdb

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
    superpose,
    tables,
)

from .attract import AttractRigidBody
from .rigidbody import RigidBody

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
    "superpose",
    "tables",
    "AttractRigidBody",
    "RigidBody",
]
