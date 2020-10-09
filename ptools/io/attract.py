
"""Read/Write Attract files."""

from ..rigidbody import RigidBody



class AttractFileParameters:
    def __init__(self, path=""):
        self.nbminim = 0
        self.minimlist = []
        self.rstk = 0.0
        self.lignames = []
        if path:
            self._init_from_file(path)

    def _init_from_file(self, path):
        """Read attract parameter file.

        Args:
            path (str): path to attract parameter file.

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
                error = "Cannot read number of minimizations from attract parameter file"
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
                error = "Cannot read minimizations from attract parameter file: "\
                        "expected {}, found {}".format(self.nbminim, i)
                raise ValueError(error) from e
            if len(tokens) < 3:
                error = "Cannot read minimization line from attract parameter file: "\
                        "expected at least 3 values, found {}".format(len(tokens))
                raise ValueError(error)
            minim = {"maxiter": int(tokens[0]),
                    "squarecutoff": float(tokens[-1]),
                    "rstk": self.rstk if tokens[-2] == '1' else 0.0}
            return minim

        # Read file ignoring comments.
        with open(path, "rt") as f:
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




def read_aminon(path):
    """Read Attract force field parameter file.

    This file is supposed to have 3 values per line:

        - index (int)
        - atom radius (float)
        - amplitude (float)
        - inull (0)

    Args:
        path (str): path to file.

    Returns:
        list[(float, float)]: radius and amplitude values for each line.
    """
    params = []
    with open(path, "rt") as f:
        for line in f:
            # Ignore empty lines and lines starting with "#".
            if line.strip() and not line.startswith("#"):
                tokens = line.split()
                radius = float(tokens[1])
                amplitude = float(tokens[2])
                params.append((radius, amplitude))
    return params


def check_ff_version_match(receptor, ligand):
    """Read force field name from receptor and ligand files and check that
    they match.

    Args:
        receptor (str): path to receptor file.
        ligand (str): path to ligand file.

    Returns:
        str: force field name.

    Raises:
        ValueError: if force field name from receptor and ligand differ.
    """
    ff_rec = read_forcefield_from_reduced(receptor)
    ff_lig = read_forcefield_from_reduced(ligand)
    if ff_rec != ff_lig:
        err = "receptor and ligand force field names do not match: '{}' != '{}'"
        err = err.format(ff_rec, ff_lig)
        raise ValueError(err)
    return ff_rec


def read_forcefield_from_reduced(path):
    """Read force field name from reduced PDB.

    Force field is read from the first line which should be formatted as
    HEADER <FORCE_FIELD_NAME>.

    Args:
        path (str): path to reduced PDB.

    Raises:
        IOError: if header cannot be extracted from first line.

    Returns:
        str: force field name in lower case.
    """
    def get_header_line():
        with open(path, "rt") as f:
            line = f.readline()
        if not line.startswith("HEADER"):
            err = (f"{path}: reduced PDB file first line must be a HEADER line"
                   "specifying the chosen force field")
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


def read_attract_parameter(path):
    return AttractFileParameters(path)


def read_translations(filename="translation.dat"):
    """Reads a translation file and returns the dictionary of translations."""
    rb = RigidBody("translation.dat")
    print(f"Read {len(rb)} translations from translation.dat")
    translations = [(atom.index, atom.coords) for atom in rb]
    return dict(translations)


def read_rotations(filename="rotation.dat"):
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
    rotdat = open("rotation.dat", "r")
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
    rotdat.close()
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
