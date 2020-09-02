
"""pyattract.io - Read/Write pyattract files."""


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
    with open(path, 'rt') as f:
        for line in f:
            # Ignore empty lines and lines starting with '#'.
            if line.strip() and not line.startswith('#'):
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
        with open(path, 'rt') as f:
            line = f.readline()
        if not line.startswith('HEADER'):
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
        except:
            error = 'Cannot read number of minimizations from attract parameter file'
            raise ValueError(error)
        return nbminim

    def read_lignames():
        def get_tokens():
            try:
                tokens = lines.pop(0).split()
            except:
                error = 'Unexpectedly reached end of attract parameter file'
                raise ValueError(error)
            return tokens
        lignames = []
        tokens = get_tokens()
        while tokens and tokens[0] == 'Lig':
            lignames.append(tokens[2])
            tokens = get_tokens()
        return lignames, tokens

    def read_rstk():
        try:
            rstk = float(tokens[3])
        except:
            error = 'Cannot read rstk from attract parameter file'
            raise ValueError(error)
        return rstk

    def read_minimization():
        try:
            tokens = lines.pop(0).split()
        except:
            error = 'Cannot read minimizations from attract parameter file: '\
                    'expected {}, found {}'.format(nbminim, i)
            raise ValueError(error)
        if len(tokens) < 3:
            error = 'Cannot read minimization line from attract parameter file: '\
                    'expected at least 3 values, found {}'.format(len(tokens))
            raise ValueError(error)
        minim = {'maxiter': int(tokens[0]),
                 'squarecutoff': float(tokens[-1]),
                 'rstk': rstk if tokens[-2] == '1' else 0.0}
        return minim

    # Read file ignoring comments.
    with open(path, 'rt') as f:
        lines = [line for line in f if not line[0] in ('#', '!')]

    # First number is the number of minimizations to perform.
    nbminim = read_number_of_minimizations()

    # Read ligand list which are all lines starting with 'Lig'.
    lignames, tokens = read_lignames()

    # Read rstk.
    rstk = read_rstk()

    # Read minimization list.
    minimlist = []
    for i in range(nbminim):
        minimlist.append(read_minimization())

    return nbminim, lignames, minimlist, rstk
