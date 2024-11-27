"""ptools.reduce.bead - Defines the ``Residue`` class."""


from typing import Any

from ..particlecollection import ParticleCollection

from .bead import Bead, BeadIdentifier
from .exceptions import (
    DuplicateBeadsError,
    MissingBeadsError,
    NoAtomsForBeadError,
    UnexpectedAtomsError,
)


class Residue:
    """Coarse-grained representation of a residue.

    Attributes:
        atoms (ParticleCollection): list of atoms that compose the residue.
        beads (list[Bead]): list of beads that compose the residue.
        reduction_parameters: list[dict[str, Any]]: reduction parameters, where each
            element is the dictionary of parameters for a bead.
    """

    name: str
    index: int
    atoms: ParticleCollection
    beads: list[Bead]
    reduction_parameters: list[dict[str, Any]]

    def __init__(
        self, atoms: ParticleCollection, bead_reduction_parameters: list[dict[str, Any]]
    ):
        """Initializes a residue."""
        self.name = atoms[0].residue_name
        self.index = atoms[0].residue_index
        self.atoms: ParticleCollection = atoms
        self.beads: list[Bead] = []
        self.reduction_parameters = bead_reduction_parameters
        self.create_beads()

    def __repr__(self) -> str:
        return f"<Residue(names={self.name!r}, {len(self.atoms)} atoms {len(self.beads)} beads)>"

    def create_beads(self):
        """Create beads from atoms and bead reduction parameters."""
        for bead_parameters in self.reduction_parameters:
            atoms = [
                atom for atom in self.atoms if atom.name in bead_parameters["atoms"]
            ]

            # Raises an error if no atoms are found for the bead.
            if not atoms:
                bead_type = bead_parameters["name"]
                bead_identifier = f"{self.name}:{self.index}:{bead_type}"
                bead_atoms = list(bead_parameters["atoms"].keys())
                raise NoAtomsForBeadError(bead_identifier, bead_atoms)

            bead = Bead(atoms, bead_parameters)
            self.beads.append(bead)

    def is_incomplete(self) -> bool:
        """Alias for ``has_missing_beads``."""
        return self.has_missing_beads()

    def has_missing_beads(self) -> bool:
        """Return ``True`` if the residue is incomplete.

        The residue is incomplete if the number of beads provided is inferior to the
        number expected from the parameters.
        """
        return len(self.beads) < len(self.reduction_parameters)

    def has_duplicate_beads(self) -> bool:
        """Return ``True`` if the residue has duplicate beads.

        The residue has duplicate beads if the number of beads provided is superior to the
        number expected from the parameters.
        """
        return len(self.beads) > len(self.reduction_parameters)

    def has_unexpected_atoms(self) -> bool:
        """Return ``True`` if the residue has atoms that are not part of any bead."""
        return len(self.atoms) > sum([len(bead.atoms) for bead in self.beads])

    def check_composition(self):
        """Check that the residue is complete, has no duplicate beads and that each
        bead's composition is correct as well."""
        for bead in self.beads:
            bead.check_composition()

        if self.has_missing_beads():
            raise MissingBeadsError(self)

        if self.has_duplicate_beads():
            raise DuplicateBeadsError(self)

        if self.has_unexpected_atoms():
            raise UnexpectedAtomsError(self)

    def find_missing_beads(self) -> set[BeadIdentifier]:
        """Return the list of missing atoms."""
        expected = {
            BeadIdentifier(
                bead_reduction_parameter["name"], bead_reduction_parameter["typeid"]
            )
            for bead_reduction_parameter in self.reduction_parameters
        }
        found = {BeadIdentifier(bead.type, bead.typeid) for bead in self.beads}
        return expected - found

    def find_duplicate_beads(self) -> set[BeadIdentifier]:
        """Return the list of duplicate atoms."""
        found = [BeadIdentifier(bead.type, bead.typeid) for bead in self.beads]
        return set([x for x in found if found.count(x) > 1])

    def find_unexpected_atoms(self) -> set[str]:
        """Return the list of unexpected atoms."""
        expected = self._expected_atoms()
        found = set([atom.name for atom in self.atoms])
        return found - expected

    def _expected_atoms(self) -> set[str]:
        """Return the list of expected atoms."""
        expected = []
        for bead in self.reduction_parameters:
            expected += list(bead["atoms"].keys())
        return set(expected)
