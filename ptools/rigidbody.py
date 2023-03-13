"""ptools.rigidbody - Defines the RigidBody class and children."""

from __future__ import annotations

# Type hinting specific imports
from collections.abc import Sequence
from typing import Optional, Type, TypeVar
from ._typing import FilePath

# Scientific libraries.
import numpy as np

# PTools.
from .atomattrs import AtomAttrs
from .io.pdb import InvalidPDBFormatError
from .io.pdb import read_pdb as io_read_pdb
from .particlecollection import ParticleCollection


RigidBodyType = TypeVar("RigidBodyType", bound="RigidBody")
AttractRigidBodyType = TypeVar("AttractRigidBodyType", bound="AttractRigidBody")


class RigidBody(ParticleCollection):
    """RigidBody is an ParticleCollection that can be initialized from a PDB file.

    It can be initialized from a file.
    """

    def __init__(self, atoms: Optional[Sequence[AtomAttrs]] = None):
        if isinstance(atoms, str):
            raise TypeError(
                "RigidBody class can not longer be instantiated from a path. "
                "Use RigidBody.from_pdb instead."
            )
        ParticleCollection.__init__(self, atoms)

    @classmethod
    def from_pdb(cls: Type[RigidBodyType], path: FilePath) -> RigidBodyType:
        atom_container = io_read_pdb(path)
        assert isinstance(atom_container, ParticleCollection)
        return cls.from_properties(atom_container.atom_properties)


class AttractRigidBody(RigidBody):
    """AttractRigidBody is a RigidBody on which one can calculate the energy.

    It has 3 additionnal arrays compared to ParticleCollection:
        - atom_categories (np.ndarray(N, )):
            1 x N shaped array for atom categories
        - atom_charges (np.ndarray(N, )):O
            1 x N shaped array for atom charges
        - atom_radii (np.ndarray(N, )):O
            1 x N shaped array for atom radii
        - atom_forces (np.ndarray(N, 3)):
            N x 3 shaped array for atom forces

    Atom categories and charges are parsed from input PDB file.
    """

    def __init__(self, atoms: Optional[Sequence[AtomAttrs]] = None):
        super().__init__(atoms)
        self._initialize_attract_properties()

    def _initialize_attract_properties(self):
        """Initializes atom categories, charges and forces from PDB extra field."""
        n_atoms = len(self)
        self.add_atom_property("category", "categories", np.zeros(n_atoms, dtype=int))
        self.add_atom_property("charge", "charges", np.zeros(n_atoms, dtype=int))
        self.add_atom_property("radius", "radii", np.zeros(n_atoms, dtype=int))
        self.add_atom_property("force", "forces", np.zeros((n_atoms, 3), dtype=int))

    def _initialize_attract_properties_from_red_extra(self):
        """Initializes atom categories, charges and forces from PDB extra field."""
        extra = self._parse_extra_from_atoms()
        self.add_atom_property(
            "category", "categories", [int(tokens[0]) - 1 for tokens in extra]
        )
        self.add_atom_property(
            "charge", "charges", [float(tokens[1]) for tokens in extra]
        )
        self.add_atom_property(
            "radius", "radii", np.zeros(len(self), dtype=float)
        )
        self.add_atom_property("force", "forces", np.zeros((len(self), 3), dtype=float))

    @classmethod
    def from_red(
        cls: Type[AttractRigidBodyType], path: FilePath
    ) -> AttractRigidBodyType:
        rigid = super().from_pdb(path)
        rigid._initialize_attract_properties_from_red_extra()
        return rigid

    @classmethod
    def from_pdb(
        cls: Type[AttractRigidBodyType], path: FilePath
    ) -> AttractRigidBodyType:
        raise NotImplementedError("Use AttractRigidBody.from_red instead.")

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
        self.forces.fill(0)

    def apply_forces(self, forces: np.ndarray):
        """Adds forces to atoms."""
        self.forces += forces
