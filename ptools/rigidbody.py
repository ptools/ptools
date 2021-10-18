"""ptools.rigidbody - Defines the RigidBody class and children."""

from __future__ import annotations
from typing import Sequence, Union

import numpy as np

from .atom import Atom, AtomCollection
from .io import read_pdb, InvalidPDBFormatError


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
        extra = self._parse_extra_from_atoms()
        self._init_atom_categories(extra)
        self._init_atom_charges(extra)

    def _parse_extra_from_atoms(self):
        """Parses extra atom field.

        Raises:
            InvalidPDBFormatError: if category or charge cannot be found or is invalid type
        """

        def assert_has_valid_category_and_charge(tokens):
            if len(tokens) < 2:
                raise InvalidPDBFormatError(
                    f"Expected atom categories and charges, found {tokens}"
                )
            assert_is_valid_category(tokens[0])
            assert_is_valid_charge(tokens[1])

        def assert_is_valid_category(token):
            """Returns True if token is an integer."""
            if not token.isdigit():
                raise InvalidPDBFormatError(
                    f"Atom category expects an int, found '{token}'"
                )

        def assert_is_valid_charge(token):
            """Returns True if token is a float."""
            try:
                float(token)
            except ValueError as error:
                raise InvalidPDBFormatError(
                    f"Atom charge expects a float, found '{token}'"
                ) from error

        extra = [atom.meta["extra"].split() for atom in self]
        for tokens in extra:
            assert_has_valid_category_and_charge(tokens)
        return extra

    def _init_atom_categories(self, extra):
        """Initializes atom category array."""
        self.atom_categories = np.array([int(tokens[0]) - 1 for tokens in extra])

    def _init_atom_charges(self, extra):
        self.atom_charges = np.array([float(tokens[1]) for tokens in extra])

    def reset_forces(self):
        """Set all atom forces to (0, 0 0)."""
        self.atom_forces.fill(0)

    def apply_forces(self, forces: np.ndarray):
        """Adds forces to atoms."""
        self.atom_forces += forces
