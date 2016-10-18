
"""pyptools.rigidbody - Defines the RigidBody class and children."""

from .atom import AtomCollection
from .io import read_pdb


class RigidBody(AtomCollection):
    """RigidBody is basically an AtomCollection that can be initialized
    from a file.

    Args:
        filename (str): path to topology file.
    """
    def __init__(self, filename):
        atoms = read_pdb(filename)
        super().__init__(atoms)
