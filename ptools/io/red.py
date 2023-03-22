"""Provides functions to read the RED file format.

The RED file format is essentially a PDB file with atoms categories and charges
in the extra field.
"""

import numpy as np

from ..particlecollection import ParticleCollection
from .._typing import FilePath
from .exceptions import InvalidREDFormatError
from .pdb import read_single_model_pdb


class _ExtraField:
    """Extra field of a PDB file.

    Provides helper functions to retrieve an atom category and charge.
    """

    def __init__(self, extra: str):
        self.extra = extra
        self.tokens = extra.split()
        self._check_number_of_tokens()

    def _check_number_of_tokens(self):
        """Checks that the extra field contains at least 2 tokens.

        If not, raises an InvalidREDFormatError.
        """
        if len(self.tokens) < 2:
            raise InvalidREDFormatError(
                f"Expected atom categories and charges, found {self.tokens!r}"
            )

    def category(self) -> int:
        """Returns the atom category."""
        try:
            category = int(self.tokens[0]) - 1
        except ValueError as error:
            raise InvalidREDFormatError(
                f"Atom category expects an int, found '{self.tokens[0]}'"
            ) from error
        return category

    def charge(self) -> float:
        """Returns the atom charge."""
        try:
            charge = float(self.tokens[1])
        except ValueError as error:
            raise InvalidREDFormatError(
                f"Atom charge expects a float, found '{self.tokens[1]}'"
            ) from error
        return charge


def read_red(path: FilePath) -> ParticleCollection:
    """Reads a RED file.

    Args:
        path: Path to the RED file.

    Returns:
        A ParticleCollection with the atoms from the RED file.

        The ParticleCollection is a regular ParticleCollection with the following
        additional properties:
          - categories: Atom categories,
          - charges: Atom charges,
    """
    atoms = read_single_model_pdb(path)
    _initialize_properties_from_extra(atoms)
    return atoms


def _initialize_properties_from_extra(atoms: ParticleCollection):
    """Initializes atom categories, charges and forces from PDB extra field."""
    extra = [_ExtraField(atom.meta["extra"]) for atom in atoms]

    natoms = len(atoms)
    atoms.add_atom_property("category", "categories", [tok.category() for tok in extra])
    atoms.add_atom_property("charge", "charges", [tok.charge() for tok in extra])
    atoms.add_atom_property("force", "forces", np.zeros((natoms, 3), dtype=float))
