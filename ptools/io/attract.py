"""Read/Write Attract files."""

import numpy as np

from ..rigidbody import RigidBody
from .._typing import FilePath


class AttractFileParameters:
    """Stores parameters from an Attract parameter file."""

    def __init__(self, path: FilePath = ""):
        self.nbminim: int = 0
        self.minimlist: list[dict[str, (int | float)]] = []
        self.rstk: float = 0.0
        self.lignames: list[str] = []
        if path:
            self._init_from_file(path)

    def _init_from_file(self, path: FilePath):
        """Read attract parameter file.

        Args:
            path (FilePath): path to attract parameter file.

        Returns:
            nbminim (int): number of minimizations to perform.
            lignames (list[str])
            minimlist (list[dict[str, (int or float)]])
            rstk (float)
        """

        def read_number_of_minimizations():
            try:
                nbminim = int(lines.pop(0).split()[0])
            except Exception as e:
                error = (
                    "Cannot read number of minimizations from attract parameter file"
                )
                raise ValueError(error) from e
            return nbminim

        def read_lignames():
            def get_tokens():
                try:
                    tokens = lines.pop(0).split()
                except Exception as e:
                    error = "Unexpectedly reached end of attract parameter file"
                    raise ValueError(error) from e
                return tokens

            lignames = []
            tokens = get_tokens()
            while tokens and tokens[0] == "Lig":
                lignames.append(tokens[2])
                tokens = get_tokens()
            return lignames, tokens

        def read_rstk():
            try:
                rstk = float(tokens[3])
            except Exception as e:
                error = "Cannot read rstk from attract parameter file"
                raise ValueError(error) from e
            return rstk

        def read_minimization():
            try:
                tokens = lines.pop(0).split()
            except Exception as e:
                error = (
                    "Cannot read minimizations from attract parameter file: "
                    f"expected {self.nbminim}, found {i}"
                )
                raise ValueError(error) from e
            if len(tokens) < 3:
                error = (
                    "Cannot read minimization line from attract parameter file: "
                    f"expected at least 3 values, found {len(tokens)}"
                )
                raise ValueError(error)
            minim = {
                "maxiter": int(tokens[0]),
                "squarecutoff": float(tokens[-1]),
                "rstk": self.rstk if tokens[-2] == "1" else 0.0,
            }
            return minim

        # Read file ignoring comments.
        with open(path, "rt", encoding="utf-8") as f:
            lines = [line for line in f if not line[0] in ("#", "!")]

        # First number is the number of minimizations to perform.
        self.nbminim = read_number_of_minimizations()

        # Read ligand list which are all lines starting with 'Lig'.
        self.lignames, tokens = read_lignames()

        # Read rstk.
        self.rstk = read_rstk()

        # Read minimization list.
        self.minimlist = []
        for i in range(self.nbminim):
            self.minimlist.append(read_minimization())


def read_aminon(path: FilePath) -> list[tuple[float, float]]:
    """Read Attract force field parameter file.

    This file is supposed to have 3 values per line:

        - index (int)
        - atom radius (float)
        - amplitude (float)
        - inull (0)

    Args:
        path (FilePath): path to file.

    Returns:
        list[(float, float)]: radius and amplitude values for each line.
    """
    params = []
    with open(path, "rt", encoding="utf-8") as f:
        for line in f:
            # Ignore empty lines and lines starting with "#".
            if line.strip() and not line.startswith("#"):
                tokens = line.split()
                radius = float(tokens[1])
                amplitude = float(tokens[2])
                params.append((radius, amplitude))
    return params


def check_ff_version_match(receptor_path: FilePath, ligand_path: FilePath) -> str:
    """Read force field name from receptor and ligand files and check that
    they match.

    Args:
        receptor_path (FilePath): path to receptor file.
        ligand_path (FilePath): path to ligand file.

    Returns:
        str: force field name.

    Raises:
        ValueError: if force field name from receptor and ligand differ.
    """
    ff_rec = read_forcefield_from_reduced(receptor_path)
    ff_lig = read_forcefield_from_reduced(ligand_path)
    if ff_rec != ff_lig:
        err = "receptor and ligand force field names do not match: '{}' != '{}'"
        err = err.format(ff_rec, ff_lig)
        raise ValueError(err)
    return ff_rec


def read_forcefield_from_reduced(path: FilePath) -> str:
    """Read force field name from reduced PDB.

    Force field is read from the first line which should be formatted as
    HEADER <FORCE_FIELD_NAME>.

    Args:
        path (FilePath): path to reduced PDB.

    Raises:
        IOError: if header cannot be extracted from first line.

    Returns:
        str: force field name in lower case.
    """

    def get_header_line():
        with open(path, "rt", encoding="utf-8") as f:
            line = f.readline()
        if not line.startswith("HEADER"):
            err = (
                f"{path}: reduced PDB file first line must be a HEADER line"
                "specifying the chosen force field"
            )
            raise IOError(err)
        return line

    def get_header_tokens():
        line = get_header_line()
        tokens = line.split()
        if len(tokens) < 2:
            err = "{}: cannot read force field name from first line '{}'"
            err = err.format(path, line)
            raise IOError(err)
        return tokens

    def get_ffname():
        ff = get_header_tokens()[1].lower()
        # if ff not in ATTRACT_FORCEFIELDS:
        #     err = "{}: invalid force field name '{}': must choose between {}"
        #     err = err.format(path, ff, ATTRACT_FORCEFIELDS)
        #     raise ValueError(err)
        return ff

    return get_ffname()


def read_attract_parameter(path: FilePath) -> AttractFileParameters:
    """Reads an Attract parameter file."""
    return AttractFileParameters(path)


def read_translations(path: FilePath = "translation.dat") -> dict[int, np.ndarray]:
    """Reads a translation file and returns the dictionary of translations."""
    rb = RigidBody.from_pdb(path)
    print(f"Read {len(rb)} translations from {path}")
    translations = [(atom.index, atom.coords) for atom in rb]
    return dict(translations)


# pylint: disable=R0914
# Not so many local variables
def read_rotations(
    path: FilePath = "rotation.dat",
) -> dict[int, tuple[float, float, float]]:
    """Returns the  dictionary of rotations read from file.

    Each rotation is a tuple (phi, theta, chi).

    The rotations in the dictionary are keyed by their file line number
    (1-based)."""
    twopi = 2.0 * 3.14159265
    nrot_per_trans = 0
    theta = []
    nphi = []
    # read theta, phi, rot data
    # nchi is number of steps to rotate about the axis joining the ligand/receptor centers
    with open(path, "rt", encoding="utf-8") as rotdat:
        line = rotdat.readline().split()
        ntheta = int(line[0])
        nchi = int(line[1])
        print(f"ntheta, nchi: {ntheta} {nchi}")
        for i in range(ntheta):
            line = rotdat.readline().split()
            theta.append(float(line[0]))
            nphi.append(int(line[1]))
            nrot_per_trans += nphi[i] * nchi
            theta[i] = twopi * theta[i] / 360.0
            print(f"{theta[i]} {nphi[i]}")

    rotations = []

    print(f"Read {ntheta} rotation lines from rotation.dat")
    print(f"{nrot_per_trans} rotations per translation")

    rotnb = 0
    for kkk in range(ntheta):
        ssii = theta[kkk]
        phii = twopi / nphi[kkk]
        for jjj in range(nphi[kkk]):
            phiii = (jjj + 1) * phii
            for iii in range(nchi):
                rotnb += 1
                chi = (iii + 1) * twopi / nchi
                rotations.append((rotnb, (phiii, ssii, chi)))

    return dict(rotations)
