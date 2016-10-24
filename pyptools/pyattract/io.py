
"""pyattract.io - Read/Write pyattract files."""


from . import PYATTRACT_FORCEFIELDS
from ..io import assert_file_exists


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
        assert_file_exists(path)
        with open(path, 'rt') as f:
            line = f.readline()
        if not line.startswith('HEADER'):
            err = '{}: reduced PDB file first line must be a HEADER line '\
                  'specifying the chosen force field ({})'
            err = err.format(path, PYATTRACT_FORCEFIELDS)
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
        if ff not in PYATTRACT_FORCEFIELDS:
            err = "{}: invalid force field name '{}': must choose between {}"
            err = err.format(path, ff, PYATTRACT_FORCEFIELDS)
            raise ValueError(err)
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
