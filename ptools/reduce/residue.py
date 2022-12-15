"""ptools.reduce.bead - Defines the ``Residue`` class."""

from dataclasses import dataclass, field, InitVar
from typing import Any

from .bead import Bead



class Residue:
    """Coarse-grained representation of a residue.

    Attributes:
        beads (list[Bead]): list of beads that compose the residue.
        residue_type (str): residue type.
        residue_index (int): residue index.
    """

    def __init__(self, residue_type: str, residue_index: int, atoms: list, bead_reduction_parameters: dict[str, Any]):
        """Initializes a residue."""
        self.beads: list[Bead] = []
        self.residue_type: str = residue_type
        self.residue_index: int = residue_index

        for bead_parameters in bead_reduction_parameters:
            atoms = [atom for atom in atoms if atom.atom_type in bead_parameters["atoms"]]



            bead = Bead(self, bead_parameters)
            self.beads.append(bead)


    def is_incomplete(self) -> bool:
        """Alias for ``has_missing_beads``."""
        return self.has_missing_beads()

    def has_missing_beads(self) -> bool:
        """Return ``True`` if the residue is incomplete.

        The residue in incomplete if the number of beads provided is inferior to the
        number expected from the parameters.
        """
        return len(self.beads) < len(self.bead_reduction_parameters)

    def has_duplicate_beads(self) -> bool:
        """Return ``True`` if the residue has duplicate beads.

        The residue has duplicate beads if the number of beads provided is superior to the
        number expected from the parameters.
        """
        return len(self.beads) > len(self.bead_reduction_parameters)

    def create_beads(self, bead_reduction_parameters: dict[str, Any]):
        """Create beads from the given parameters.

        Args:
            bead_reduction_parameters: reduction parameters, where keys are parameters names
                and values are parameters values.
        """
        self.bead_reduction_parameters = bead_reduction_parameters
        self.beads = []
        for bead_reduction_parameters in self.bead_reduction_parameters.values():
            beads = self._create_beads(bead_reduction_parameters)
            self.beads.extend(beads)

    def _create_beads(self, bead_reduction_parameters: dict[str, Any]):
        """Create beads from the given parameters.

        Args:
            bead_reduction_parameters: reduction parameters, where keys are parameters names
                and values are parameters values.

        Returns:
            list of beads.
        """
        beads = []
        for atom_reduction_parameters in bead_reduction_parameters["atoms"].values():
            atoms = self._get_atoms(atom_reduction_parameters)
            if len(atoms) == 0:
                continue
            bead = Bead(atoms, self.resid_type, self.resid_id, bead_reduction_parameters)
            beads.append(bead)
        return beads


