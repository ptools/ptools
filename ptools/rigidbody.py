"""ptools.rigidbody - Defines the RigidBody class and children."""

from __future__ import annotations

from abc import ABC, abstractmethod
import os
from typing import Sequence

import numpy as np

from .atom import Atom, AtomCollection
from .io.pdb import InvalidPDBFormatError
from .io.pdb import read_pdb as io_read_pdb


class _RigidBodyBase(ABC):
    """Base Base class for RigidBody subclasses.

    Ensures that child classes implement the `read_pdb` method.
    """

    @abstractmethod
    def read_pdb(self, path: str | bytes | os.PathLike):
        """Initializes internal components from data read in PDB file."""


class RigidBodyBase(_RigidBodyBase, AtomCollection):
    """Base class for RigidBody subclasses.

    Implements RigidBodyBase.from_pdb which returns a new instance initialized
    from a PDB file, read with `cls.read_pdb`.
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
        """Post initialization method.

        Overriden if needed in child classes.
        """

    @classmethod
    def from_pdb(cls, path: str | bytes | os.PathLike):
        rigid = cls()
        rigid.read_pdb(path)
        return rigid


class RigidBody(RigidBodyBase):
    """RigidBody is basically an AtomCollection that can be initialized
    from a file.
    """

    def read_pdb(self, path: str | bytes | os.PathLike):
        atoms = io_read_pdb(path)
        super().__init__(atoms)


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

    def __post_init__(self):
        """Read each atom 'extra' meta-data attribute to find atom charge
        and category.

        Raises:
            IOError: if atom meta['extra'] attribute does not contain 2 entries.
        """
        N = len(self)
        self.atom_categories = np.zeros(N, dtype=int)
        self.atom_charges = np.zeros(N, dtype=float)
        self.atom_forces = np.zeros((N, 3), dtype=float)

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
