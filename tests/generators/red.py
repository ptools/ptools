"""Helper classes and functions for RED format I/O testing."""

from dataclasses import dataclass
import tempfile

from ptools.tables import atomic_radii

from .io import generate_tmp_file


@dataclass
class RedFileBuilder:
    """Creates a standard RED file."""

    has_categories: bool = True
    has_charges: bool = True
    has_invalid_charges: bool = False

    @staticmethod
    def base_atom_list() -> list[str]:
        """Returns basic atoms (no categories, no charges)."""
        return [
            "ATOM      1  CA  TYR     1       0.880   6.552  -1.114",
            "ATOM      2  CSE TYR     1      -0.313   5.011  -0.781",
            "ATOM      3  CSE TYR     1      -0.074   2.173  -0.857",
            "ATOM      4  CA  SER     2       3.501   7.012   1.646",
            "ATOM      5  CSE SER     2       3.640   8.268   1.516",
            "ATOM      6  CA  SER     3       1.143   5.686   4.387",
            "ATOM      7  CSE SER     3       0.555   6.008   5.467",
            "ATOM      8  CA  TYR     4      -1.411   2.949   3.745",
            "ATOM      9  CSE TYR     4      -1.675   0.958   3.911",
            "ATOM     10  CSE TYR     4      -3.525  -1.079   3.005",
        ]

    @property
    def _natoms(self) -> int:
        """Number of atoms in the RED file."""
        return len(self.base_atom_list())

    def atoms(self) -> str:
        """Atoms in RED format."""
        atoms = self.base_atom_list()

        if self.has_categories:
            categories = self.categories()
            for i, atom in enumerate(atoms):
                atoms[i] = atom + f"{categories[i]:5d}"

        if self.has_charges:
            if self.has_invalid_charges:
                charges = self.invalid_charges()
            else:
                charges = [f"{charge:8.3f}" for charge in self.charges()]
            for i, atom in enumerate(atoms):
                atoms[i] = atom + f"{charges[i]:>8s}"

        return "\n".join(atoms)

    @classmethod
    def categories(cls) -> list[int]:
        """Atom categories as presented in default atoms."""
        return [1, 27, 28, 1, 23, 1, 23, 1, 27, 28]

    @classmethod
    def charges(cls) -> list[float]:
        """Atom charges as presented in default atoms."""
        return [0.0] * len(cls.base_atom_list())

    @classmethod
    def radii(cls) -> list[float]:
        """Atom radii as presented in default atoms."""
        return [atomic_radii["C"]] * len(cls.base_atom_list())

    @classmethod
    def invalid_charges(cls) -> list[str]:
        """Atom invalid charges (strings instead of floats)."""
        return ["FOO"] * len(cls.base_atom_list())

    @property
    def content(self) -> str:
        """Returns the content of the RED file."""
        return self.atoms()


def generate_red_file(
    has_categories: bool = True, has_charges: bool = True, invalid_charges: bool = False
) -> tempfile.NamedTemporaryFile:
    """Creates a temporary file (RED format) that contains 10 atoms."""
    if invalid_charges:
        return generate_tmp_file(
            content=RedFileBuilder(has_invalid_charges=True).content
        )
    return generate_tmp_file(
        content=RedFileBuilder(has_categories, has_charges).content
    )
