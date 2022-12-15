"""ptools.reduce.exceptions - Exceptions for ``ptools.reduce``."""

import collections


class ResidueReductionError(Exception):
    """Base class for exceptions raised when an error occurs during
    a residue reduction (transformation from all-atom to coarse grain).

    Attrs:
        default_message_fmt (str): default message format.

    Args:
        residue_name (str): residue name.
        residue_id (int): residue identifier.
        kwargs dict[str]->value: any argument that will be used to format
            message.
    """

    default_message_fmt = "in residue {self.residue_name}:{self.residue_id}"

    def __init__(self, residue_name, residue_id, message="", **kwargs):
        self.residue_name = residue_name
        self.residue_id = residue_id

        # Sets attributes from kwargs.
        for attr, value in kwargs.items():
            if attr not in ("residue_name", "residue_id", "message"):
                setattr(self, attr, value)

        message = message or self.default_message_fmt
        super().__init__(message.format(self=self))

    def report(self):
        return f"{type(self).__name__}: {self.message}"


class NoResidueReductionRulesFoundError(ResidueReductionError):
    """Raised when no residue reduction rules are found for a given residue."""

    default_message_fmt = (
        ResidueReductionError.default_message_fmt
        + ": no rule found for this residue reduction"
    )


class BeadCreationError(ResidueReductionError):
    """Base class for exceptions raised when an error is encountered when reducing
    an atomistic topology to a coarse grain topology.

    Attrs:
        bead (ptools.commands.reduce_cmd.Bead): bead where error occured.
        expected_atoms (list[str]): list of expected atom names.
        found_atoms (list[str]): list of actually found atom names.
    """

    default_message_fmt = (
        ResidueReductionError.default_message_fmt + " in bead '{self.bead.name}'"
    )

    def __init__(self, bead):
        self.bead = bead
        self.expected_atoms = sorted(self.bead.atom_reduction_parameters.keys())
        self.found_atoms = sorted(atom.atom_type for atom in self.bead.atoms)
        super().__init__(bead.resid_type, bead.resid_id, bead_name=self.bead.name)


class IncompleteBeadError(BeadCreationError):
    """Raised when an atom is missing when constructing a bead."""

    default_message_fmt = (
        BeadCreationError.default_message_fmt
        + "\n    expected atoms: {self.expected_atoms}"
        "\n    found atoms: {self.found_atoms}"
        "\n    missing atoms: {self.missing_atoms}"
    )

    @property
    def missing_atoms(self):
        """Return the list of missing atom names."""
        return list(set(self.expected_atoms) - set(self.found_atoms))


class DuplicateAtomInBeadError(BeadCreationError):
    """Raised when an atom is found several times when constructing a bead."""

    default_message_fmt = (
        BeadCreationError.default_message_fmt
        + "\n    expected atoms: {self.expected_atoms}"
        "\n    found atoms: {self.found_atoms}"
        "\n    duplicate atoms: {self.duplicate_atoms}"
    )

    @property
    def duplicate_atoms(self) -> list[str]:
        """Return the list of duplicate atom names."""
        counter = collections.Counter(self.found_atoms)
        return [name for name, count in counter.items() if count > 1]


class IgnoredAtomsInReducedResidueError(ResidueReductionError):
    """Raised when an atom from all-atom model has not been used in the coarse grain model."""

    default_message_fmt = (
        ResidueReductionError.default_message_fmt
        + " some atoms were unused during coarse grain "
        "modeling: {self.unused_atoms}"
    )

    def __init__(self, residue):
        # Names of atoms sent to coarse grain residue constructor.
        self.atom_names = [atom.atom_type for atom in residue.all_atoms]

        # Names of atoms actually used in the beads.
        self.bead_atom_names = [
            atom.atom_type for bead in residue.beads for atom in bead.atoms
        ]

        super().__init__(residue.resname, residue.resid)

    @property
    def unused_atoms(self) -> list[str]:
        """Returns the list of unused atom names."""
        diff = set(self.atom_names) ^ set(self.bead_atom_names)
        return list(diff)
