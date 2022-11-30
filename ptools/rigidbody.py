"""ptools.rigidbody - Defines the RigidBody class and children."""

from __future__ import annotations

# Type hinting specific imports
from collections.abc import Sequence
from typing import Type, TypeVar
from ._typing import FilePath

# Scientific libraries.
import numpy as np

# PTools.
from .atomattrs import AtomAttrs
from .atomcollection import AtomCollection
from .io.pdb import InvalidPDBFormatError
from .io.pdb import read_pdb as io_read_pdb


RigidBodyType = TypeVar("RigidBodyType", bound="RigidBody")
AttractRigidBodyType = TypeVar("AttractRigidBodyType", bound="AttractRigidBody")


class RigidBody(AtomCollection):
    """RigidBody is an AtomCollection that can be initialized from a PDB file.

    It has 3 additionnal arrays comparaed to AtomCollection:
        - atom_categories (np.ndarray(N, )):
            1 x N shaped array for atom categories
        - atom_charges (np.ndarray(N, )):O
            1 x N shaped array for atom charges
        - atom_forces (np.ndarray(N, 3)):
            N x 3 shaped array for atom forces

    Also, it can be initialized from a file.
    """

    def __init__(self, atoms: Sequence[AtomAttrs] = None):
        if isinstance(atoms, str):
            raise TypeError(
                "RigidBody class can not longer be instantiated from a path. "
                "Use RigidBody.from_pdb instead."
            )
        AtomCollection.__init__(self, atoms)

    @classmethod
    def from_pdb(cls: Type[RigidBodyType], path: FilePath) -> RigidBodyType:
        atoms = io_read_pdb(path)
        assert isinstance(atoms, AtomCollection)
        return cls(atoms)


class AttractRigidBody(RigidBody):
    """AttractRigidBody is a RigidBody on which one can calculate the energy.

    Atom categories and charges are parsed from input PDB file."""

    def __init__(self, atoms: Sequence[AtomAttrs] = None):
        super().__init__(atoms)
        N = len(self)
        self.atom_categories = np.zeros(N, dtype=int)
        self.atom_charges = np.zeros(N, dtype=float)
        self.atom_forces = np.zeros((N, 3), dtype=float)

    @classmethod
    def from_pdb(
        cls: Type[AttractRigidBodyType], path: FilePath
    ) -> AttractRigidBodyType:
        rigid = super().from_pdb(path)
        rigid._init_categories_and_charges_from_pdb_extra()
        return rigid

    def _init_categories_and_charges_from_pdb_extra(self):
        """Initializes atom and charges arrays from data read from PDB."""
        extra = self._parse_extra_from_atoms()
        self.atom_categories = np.array([int(tokens[0]) - 1 for tokens in extra])
        self.atom_charges = np.array([float(tokens[1]) for tokens in extra])
        self.atom_forces = np.zeros((len(self), 3), dtype=float)

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

    def reset_forces(self):
        """Set all atom forces to (0, 0 0)."""
        self.atom_forces.fill(0)

    def apply_forces(self, forces: np.ndarray):
        """Adds forces to atoms."""
        self.atom_forces += forces
