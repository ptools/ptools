"""Helper classes and functions for RED format I/O testing."""

import tempfile
from .temporaryfile import mk_tmp_file

class TestREDBuilder:
    """Creates a standard RED file."""

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

    @classmethod
    def atom(cls, has_categories: bool = True, has_charges: bool = True) -> str:
        """Atoms in RED format."""
        atoms = cls.base_atom_list()

        if has_categories:
            categories = cls.categories()
            for i, atom in enumerate(atoms):
                atoms[i] = atom + f"{categories[i]:5d}"

        if has_charges:
            charges = cls.charges()
            for i, atom in enumerate(atoms):
                atoms[i] = atom + f"{charges[i]:8.3f}"

        return "\n".join(atoms)

    @classmethod
    def categories(cls) -> list[int]:
        """Atom categories as presented in default atoms."""
        return [1, 27, 28, 1, 23, 1, 23, 1, 27, 28]

    @classmethod
    def charges(cls) -> list[float]:
        """Atom charges as presented in default atoms."""
        return [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]

    @classmethod
    def redfile(cls, has_categories: bool = True, has_charges: bool = True) -> str:
        return cls.atom(has_categories, has_charges)

    @classmethod
    def redfile_invalid_charges(cls) -> str:
        atoms = cls.redfile(has_charges=False).splitlines()
        for i, atom in enumerate(atoms):
            atoms[i] = atom + f"{'FOO':>8s}"
        return "\n".join(atoms)


def mk_red(has_categories: bool = True, has_charges: bool = True) -> tempfile.NamedTemporaryFile:
    """Creates a temporary file (RED format) that contains 10 atoms."""
    return mk_tmp_file(content=TestREDBuilder.redfile(has_categories, has_charges))


def mk_red_invalid_charges() -> tempfile.NamedTemporaryFile:
    """Creates a temporary file (RED format) that contains 10 atoms and dummy charges."""
    return mk_tmp_file(content=TestREDBuilder.redfile_invalid_charges())
