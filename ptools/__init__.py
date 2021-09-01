"""Top-level package for PTools"""

__author__ = """Benoist LAURENT"""
__email__ = "benoistlaurent@gmail.com"
__version__ = "0.1.0"


from . import (
    atom,
    attract,
    forcefield,
    heligeom,
    io,
    pairlist,
    rigidbody,
    spatial,
    superpose,
    tables,
)

from .rigidbody import AttractRigidBody, RigidBody
