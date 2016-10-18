

import os
import tempfile

import numpy


TEST_PDB = os.path.join(os.path.dirname(__file__), 'data', 'test_10atoms.pdb')
TEST_RED = os.path.join(os.path.dirname(__file__), 'data', 'test_10atoms.red')
TEST_PDB_ATOM_NAMES = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ', 'H1']

assert_array_almost_equal = numpy.testing.assert_array_almost_equal


def random_filename():
    """Return a random file name."""
    tmpfile = tempfile.NamedTemporaryFile()
    tmpfile.close()
    return tmpfile.name


def mk_tmp_file(content='', mode='wt', **kwargs):
    """Create a temporary empty file.
    Returns:
        tempfile.NamedTemporaryFile
    """
    tmpfile = tempfile.NamedTemporaryFile(mode, **kwargs)
    tmpfile.write(content)
    tmpfile.flush()
    tmpfile.seek(0)
    return tmpfile


def mk_empty_file(**kwargs):
    return mk_tmp_file(content='', **kwargs)
