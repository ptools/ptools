"""ptools.rigidbody - Defines the RigidBody class and children."""

from __future__ import annotations
from typing import Sequence, Union

import numpy as np

from .atom import Atom, AtomCollection
from .io import read_pdb


class RigidBody(AtomCollection):
    """RigidBody is basically an AtomCollection that can be initialized
    from a file.
    """

    def __init__(self, atoms_or_path: Union[Sequence[Atom], str] = None):
        if isinstance(atoms_or_path, str):
            atoms_or_path = read_pdb(atoms_or_path)

        # At this point, atoms_or_path is list[Atom] | None
        super().__init__(atoms_or_path)


class AttractRigidBody(RigidBody):
    """AttractRigidBody is a RigidBody with several force field features.

    Attributes:
        with `N` being the number of atoms

        atom_categories (np.ndarray(N, )):
            1 x N shaped array for atom categories
        atom_charges (np.ndarray(N, )):
            1 x N shaped array for atom charges
        atom_forces (np.ndarray(N, 3)):
            N x 3 shaped array for atom forces
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        N = len(self)
        self.atom_categories = np.zeros(N, dtype=int)
        self.atom_charges = np.zeros(N, dtype=float)
        self.atom_forces = np.zeros((N, 3), dtype=float)
        self._init()

    def _init(self):
        """Read each atom 'extra' meta-data attribute to find atom charge
        and category.

        Raises:
            IOError: if atom meta['extra'] attribute does not contain 2 entries.
        """

        def init_category():
            try:
                self.atom_categories = np.array(
                    [int(tokens[0]) - 1 for tokens in extra]
                )
            except Exception as e:
                err = f"cannot initialize atom category: {e}"
                raise IOError(err) from e

        def init_charges():
            try:
                self.atom_charges = np.array([float(tokens[1]) for tokens in extra])
            except Exception as e:
                err = f"cannot initialize atom charges: {e}"
                raise IOError(err) from e

        extra = [atom.meta["extra"].split() for atom in self]
        init_category()
        init_charges()

    def reset_forces(self):
        """Set all atom forces to (0, 0 0)."""
        self.atom_forces.fill(0)

    def apply_forces(self, forces: np.ndarray):
        """Adds forces to atoms."""
        self.atom_forces += forces
