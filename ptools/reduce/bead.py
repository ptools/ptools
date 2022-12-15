"""ptools.reduce.bead - Defines the ``Bead`` class."""

from typing import Any, Protocol

from .exceptions import IncompleteBeadError, DuplicateAtomInBeadError



class Bead:
    """"A bead is a coarse-grained representation of several atoms.

    Attributes:
        atoms: list of atoms that compose the bead.
        residue_type: residue type.
        residue_index: residue index.
        type: bead type.
        charge: bead charge.
        typeid: bead type id (present in PDB format in 'extra' field).
    """

    atoms: list
    residue_type: str
    residue_index: int
    type: str
    charge: float
    typeid: int
    chain: str
    atom_reduction_parameters: dict[str, Any]

    def __init__(self, atoms, residue_type: str, residue_index: int, bead_reduction_parameters: dict[str, Any]):
        """Initialize a bead.

        Args:
            atoms: list of atoms that compose the bead.
            bead_reduction_parameters: reduction parameters, where keys are parameters names
                and values are parameters values.
        """
        self.atoms = atoms
        self.residue_type = residue_type
        self.residue_index = residue_index
        self.type = bead_reduction_parameters.get("name", "X")
        self.charge = bead_reduction_parameters.get("charge", 0.0)
        self.typeid = bead_reduction_parameters.get("typeid", 0)
        self.atom_reduction_parameters = bead_reduction_parameters["atoms"]

        # Chain is set from the first atom of the bead.
        self.chain = atoms[0].chain

    def is_incomplete(self) -> bool:
        """Alias for ``has_missing_atoms``."""
        return self.has_missing_atoms()

    def has_missing_atoms(self) -> bool:
        """Return ``True`` if the bead is incomplete.

        The bead in incomplete if the number of atoms provided is inferior to the
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
        if self.has_missing_atoms():
            raise IncompleteBeadError(self)
        if self.has_duplicate_atoms():
            raise DuplicateAtomInBeadError(self)


