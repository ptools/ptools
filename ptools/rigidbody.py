"""ptools.rigidbody - Defines the RigidBody class and children."""

from __future__ import annotations

# Python core libraries.
import os

# Scientific libraries.
import numpy as np

# PTools.
from .atom import AtomCollection
from .io.pdb import InvalidPDBFormatError, FromPDB
from .io.pdb import read_pdb as io_read_pdb


class RigidBody(FromPDB, AtomCollection):
    """RigidBody is an AtomCollection on which one can calculate the energy.

    It has 3 additionnal arrays comparaed to AtomCollection:
        - atom_categories (np.ndarray(N, )):
            1 x N shaped array for atom categories
        - atom_charges (np.ndarray(N, )):
            1 x N shaped array for atom charges
        - atom_forces (np.ndarray(N, 3)):
            N x 3 shaped array for atom forces

    Also, it can be initialized from a file.
    """

    def __init__(self, *args, **kwargs):
        if args and isinstance(args[0], str):
            raise TypeError(
                "RigidBody class can not longer be instantiated from a path. "
                "Use RigidBody.from_pdb instead."
            )

        AtomCollection.__init__(self, *args, **kwargs)
        self.__post_init__()

    def __post_init__(self):
        """Post-initialization.

        Necessary to initialize a new instance using class.from_pdb.
        """
        N = len(self)
        self.atom_categories = np.zeros(N, dtype=int)
        self.atom_charges = np.zeros(N, dtype=float)
        self.atom_forces = np.zeros((N, 3), dtype=float)

    def reset_forces(self):
        """Set all atom forces to (0, 0 0)."""
        self.atom_forces.fill(0)

    def apply_forces(self, forces: np.ndarray):
        """Adds forces to atoms."""
        self.atom_forces += forces

    def init_from_pdb(self, path: str | bytes | os.PathLike):
        """Initialization from PDB file."""
        atoms = io_read_pdb(path)
        AtomCollection.__init__(self, atoms)
        self.__post_init__()


class AttractRigidBody(RigidBody):
    """AttractRigidBody is a RigidBody with custom initialization from PDB file."""

    def __post_init__(self):
        """Post-initialization.

        Necessary to initialize a new instance using class.from_pdb.
        """
        super().__post_init__()
        self._init_categories_and_charges_from_pdb_extra()

    def _init_categories_and_charges_from_pdb_extra(self):
        """Initializez atom and charges arrays from data read from PDB."""
        extra = self._parse_extra_from_atoms()
        self.atom_categories = np.array([int(tokens[0]) - 1 for tokens in extra])
        self.atom_charges = np.array([float(tokens[1]) for tokens in extra])

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
