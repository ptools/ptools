"""ptools.reduce.bead - Defines the ``Bead`` class."""

from dataclasses import dataclass
from typing import Any

import numpy as np

from .exceptions import IncompleteBeadError, DuplicateAtomsError
from ..linalg import center_of_mass
from ..particlecollection import ParticleCollection


@dataclass(frozen=True)
class BeadIdentifier:
    """Identifier for a bead.

    Attributes:
        type (str): type (name) of the bead.
        typeid (int): numeric identifier
    """

    type: str
    typeid: int

    def __hash__(self) -> int:
        return hash((self.type, self.typeid))


@dataclass
class BeadDump:
    """A bead with no reference to original atoms.

    Handy for conversion to ParticleCollection.
    """
    index: int
    name: str
    residue_name: str
    residue_index: int
    type: str
    typeid: int
    charge: float
    chain: str
    coordinates: np.ndarray


class Bead:
    """ "A bead is a coarse-grained representation of several atoms.

    Attributes:
        atoms: list of atoms that compose the bead.
        index: bead index.
        residue_name: residue type.
        residue_index: residue index.
        type: bead type.
        typeid: bead type id (present in PDB format in 'extra' field).
        charge: bead charge.
        chain: chain identifier.
    """

    atoms: ParticleCollection
    index: int
    residue_name: str
    residue_index: int
    type: str
    typeid: int
    charge: float
    chain: str
    atom_reduction_parameters: dict[str, Any]

    def __init__(self, atoms, bead_reduction_parameters: dict[str, Any]):
        """Initialize a bead.

        Args:
            atoms: list of atoms that compose the bead.
            bead_reduction_parameters: reduction parameters, where keys are parameters names
                and values are parameters values.
        """
        if len(atoms) == 0:
            raise ValueError("Bead must have at least one atom.")

        self.atoms = atoms
        self.index = 0
        self.chain = atoms[0].chain
        self.residue_name = self.atoms[0].residue_name
        self.residue_index = self.atoms[0].residue_index
        self.type = bead_reduction_parameters.get("name", "X")
        self.typeid = bead_reduction_parameters.get("typeid", 0)
        self.charge = bead_reduction_parameters.get("charge", 0.0)
        self.atom_reduction_parameters = bead_reduction_parameters["atoms"]

    def dump(self) -> BeadDump:
        """Return a bead dump."""
        return BeadDump(
            index=self.index,
            name=self.name,
            residue_name=self.residue_name,
            residue_index=self.residue_index,
            type=self.type,
            typeid=self.typeid,
            charge=self.charge,
            chain=self.chain,
            coordinates=self.coordinates,
        )

    @property
    def coordinates(self) -> np.ndarray:
        """Returns the coordinates of the bead.

        The coordinates of the bead are the center of mass of the bead atoms.
        """
        coordinates = [atom.coordinates for atom in self.atoms]
        weights = [
            self.atom_reduction_parameters[atom.name].get("weight", 1.0)
            for atom in self.atoms
        ]
        return center_of_mass(coordinates, weights)

    @property
    def name(self) -> str:
        """Alias for ``type``."""
        return self.type

    def is_incomplete(self) -> bool:
        """Alias for ``has_missing_atoms``."""
        return self.has_missing_atoms()

    def has_missing_atoms(self) -> bool:
        """Return ``True`` if the bead is incomplete.

        The bead is incomplete if the number of atoms provided is inferior to the
        number expected from the parameters.
        """
        return len(self.atoms) < len(self.atom_reduction_parameters)

    def has_duplicate_atoms(self) -> bool:
        """Return ``True`` if the bead has duplicate atoms.

        The bead has duplicate atoms if the number of atoms provided is superior to the
        number expected from the parameters.
        """
        return len(self.atoms) > len(self.atom_reduction_parameters)

    def check_composition(self):
        """Check that the bead is complete and has no duplicate atoms."""
        bead_identifier = f"{self.residue_name}:{self.residue_index}:{self.type}"
        if self.has_missing_atoms():
            raise IncompleteBeadError(bead_identifier, self.find_missing_atoms())
        if self.has_duplicate_atoms():
            raise DuplicateAtomsError(bead_identifier, self.find_duplicate_atoms())

    def find_missing_atoms(self) -> list[str]:
        """Find the missing atoms in the bead.

        Returns:
            list of missing atom names.
        """
        expected_atoms = list(self.atom_reduction_parameters.keys())
        found = [atom.name for atom in self.atoms]
        missing = [atom for atom in expected_atoms if atom not in found]
        return missing

    def find_duplicate_atoms(self) -> list[str]:
        """Find the duplicate atoms in the bead.

        Returns:
            list of duplicate atom names.
        """
        found = [atom.name for atom in self.atoms]
        duplicate = [atom for atom in found if found.count(atom) > 1]
        return list(set(duplicate))
