
"""pyptools.rigidbody - Defines the RigidBody class and children."""

import copy

import numpy

from .atom import AtomCollection
from .io import read_pdb


class RigidBody(AtomCollection):
    """RigidBody is basically an AtomCollection that can be initialized
    from a file.

    Args:
        arg (str, RigidBody): path to topology file or parent RigidBody.
    """
    def __init__(self, arg):
        if isinstance(arg, str):
            atoms = read_pdb(arg)
        elif isinstance(arg, AtomCollection):
            atoms = copy.deepcopy(arg.atoms)
        else:
            err = 'RigidBody can only be initialized from file name or '\
                  'parent RigidBody'
            raise TypeError(err)
        super().__init__(atoms)


class AttractRigidBody(RigidBody):
    """AttractRigidBody is a RigidBody with several force field features.

    Attributes:
        with `N` being the number of atoms

        atom_categories (numpy.ndarray(N, )):
            1 x N shaped array for atom categories
        atom_charges (numpy.ndarray(N, )):
            1 x N shaped array for atom charges
        atom_forces (numpy.ndarray(N, 3)):
            N x 3 shaped array for atom forces
    """

    def __init__(self, filename):
        super().__init__(filename)
        N = len(self)
        self.atom_categories = numpy.zeros(N, dtype=int)
        self.atom_charges = numpy.zeros(N, dtype=float)
        self.atom_forces = numpy.zeros((N, 3), dtype=float)
        self._init()

    def _init(self):
        """Read each atom 'extra' meta-data attribute to find atom charge
        and category.

        Raises:
            IOError: if atom meta['extra'] attribute does not contain 2 entries.
        """
        def init_category():
            try:
                self.atom_categories = numpy.array([int(tokens[0]) - 1 for tokens in extra])
            except Exception as e:
                err = 'cannot initialize atom category: {}'.format(e)
                raise IOError(err)

        def init_charges():
            try:
                self.atom_charges = numpy.array([float(tokens[1]) for tokens in extra])
            except Exception as e:
                err = 'cannot initialize atom charges: {}'.format(e)
                raise IOError(err)

        extra = [atom.meta['extra'].split() for atom in self.atoms]
        init_category()
        init_charges()

    def reset_forces(self):
        """Set all atom forces to (0, 0 0)."""
        self.atom_forces.fill(0)

