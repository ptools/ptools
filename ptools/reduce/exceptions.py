"""ptools.reduce.exceptions - Exceptions for ``ptools.reduce``."""

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from .residue import Residue


class NoReductionRulesError(ValueError):
    """Raised when no reduction rules are found for a residue.

    Args:
        residue_name (str): residue name.
        residue_id (int): residue identifier.
    """

    def __init__(self, residue_name, residue_id):
        message = f"No reduction rules found for residue {residue_name}:{residue_id}."
        super().__init__(message)


class MissingBeadsError(ValueError):
    """Raised when a residue has missing beads.

    Args:
        residue_name (str): residue name.
        residue_id (int): residue identifier.
        missing_beads (list[BeadIdentifier]): list of missing bead identifiers.
    """

    def __init__(self, residue: "Residue"):
        message = (
            f"Residue {residue.name}:{residue.index} is incomplete. "
            f"Missing beads: {residue.find_missing_beads()}."
        )
        super().__init__(message)


class DuplicateBeadsError(ValueError):
    """Raised when a residue has duplicate beads.

    Args:
        residue (Residue): residue that has duplicate beads.
    """

    def __init__(self, residue: "Residue"):
        message = (
            f"Residue {residue.name}:{residue.index} has duplicate beads. "
            f"Duplicate beads: {residue.find_duplicate_beads()}."
        )
        super().__init__(message)


class IncompleteBeadError(ValueError):
    """Exception raised when a bead is incomplete."""

    def __init__(self, bead_identifier: str, missing_atoms: list[str]):
        """Initialize the exception.

        Args:
            bead_identifier (str): bead identifier (<residue_name>:<residue_index>:<bead_name>)
            missing_atoms (list[str]): list of missing atom names.
        """
        super().__init__(
            f"Bead {bead_identifier} is incomplete, missing atoms: {missing_atoms}"
        )


class NoAtomsForBeadError(ValueError):
    """Exception raised when no atoms are found for a bead."""

    def __init__(self, bead_identifier: str, bead_atoms: list[str]):
        """Initialize the exception.

        Args:
            bead_identifier (str): bead identifier (<residue_name>:<residue_index>:<bead_name>)
            bead_atoms (list[str]): list of bead atom names.
        """
        super().__init__(
            f"Bead {bead_identifier}: no atoms found, expected: {bead_atoms}"
        )


class DuplicateAtomsError(ValueError):
    """Exception raised when a bead has duplicate atoms."""

    def __init__(self, bead_identifier: str, duplicate_atoms: list[str]):
        """Initialize the exception.

        Args:
            bead_identifier (str): bead identifier (<residue_name>:<residue_index>:<bead_name>)
            duplicate_atoms (list[str]): list of duplicate atom names.
        """
        super().__init__(
            f"Bead {bead_identifier} has duplicate atoms: {duplicate_atoms}"
        )


class UnexpectedAtomsError(Exception):
    """Exception raised when the residue has atoms that are not part of any bead."""

    def __init__(self, residue: Residue):
        """Initializes the exception."""
        self.unexpected_atoms = residue.find_unexpected_atoms()
        self.message = (
            f"Residue {residue.name}:{residue.index} has atoms that are not part of "
            f"any bead: {self.unexpected_atoms}"
        )
        super().__init__(self.message)
