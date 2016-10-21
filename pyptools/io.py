
"""linalg - Defines functions to read/write files."""

import os

from .atom import BaseAtom, AtomCollection


def read_pdb(filename):
    """Read a Protein Data Bank file.

    Args:
        filename (str): path to file.

    Returns:
        AtomCollection: collection of Atoms
    """
    def get_header(line):
        """Returns PDB line header i.e. first 6 characters."""
        return line[:6]

    def is_atom_line(line):
        """Return True is the line header is 'ATOM  '."""
        return get_header(line) == 'ATOM  '

    def read_atom_line(line):
        """Return an `atom.BaseAtom` instance initialized with data read from
        line."""
        index = int(line[6:11])
        name = line[12:16].strip().upper()
        resname = line[17:20].strip().upper()
        chain = line[21].strip()
        resid = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coords = (x, y, z)
        extra = line[54:].strip()

        atom = BaseAtom(index=index, name=name,
                        resname=resname, resid=resid,
                        chain=chain,
                        coords=coords,
                        meta={'extra': extra})
        return atom

    atoms = []
    assert_file_exists(filename)
    with open(filename, 'rt') as f:
        for line in f:
            if is_atom_line(line):
                atoms.append(read_atom_line(line))
    return AtomCollection(atoms)


def assert_file_exists(path):
    """Check input file exists, raises an exception if it does not.

    Args:
        path (str): path to file.

    Raises:
        FileNotFound: if file does not exist.
    """
    if not os.path.exists(path):
        err = "No such file: '{}'".format(path)
        raise FileNotFoundError(err)
