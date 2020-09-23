
"""test_io - Tests for `ptools.io` module."""

import unittest

import ptools.io.pdb as pdb

from . import TEST_PDB, TEST_PDB_3MODELS, TEST_PDB_ATOM_NAMES
from .testing.io import random_filename


class TestIO(unittest.TestCase):

    def test_read_pdb_does_not_exist(self):
        filename = random_filename()
        with self.assertRaises(FileNotFoundError):
            pdb.read_pdb(filename)

    def test_is_atom_line(self):
        self.assertTrue(pdb.is_atom_line("ATOM  <rest_of_the_line>"))
        self.assertTrue(pdb.is_atom_line("HETATM<rest_of_the_line>"))
        self.assertFalse(pdb.is_atom_line("ATOMI"))
        self.assertFalse(pdb.is_atom_line("ANISOU"))

    def test_read_pdb(self):
        atoms = pdb.read_pdb(TEST_PDB)
        self.assertEqual(len(atoms), 10)
        for i in range(10):
            atom = atoms[i]
            self.assertEqual(atom.index, i + 1)
            self.assertEqual(atom.name, TEST_PDB_ATOM_NAMES[i])
            self.assertEqual(atom.resid, 1)
            self.assertEqual(atom.resname, 'LYS')
            self.assertEqual(atom.chain, 'A')
            self.assertAlmostEqual(atom.coords[0], i + 1)
            self.assertAlmostEqual(atom.coords[1], 10 + i + 1)
            self.assertAlmostEqual(atom.coords[2], 20 + i + 1)

    def test_read_pdb_multiple_models(self):
        models = pdb.read_pdb(TEST_PDB_3MODELS)
        self.assertEqual(len(models), 3)
        for atoms in models:
            self.assertTrue(len(atoms), 10)



# ---------- coverage: platform darwin, python 3.8.1-final-0 -----------
# Name                    Stmts   Miss  Cover   Missing
# -----------------------------------------------------
# ptools/__init__.py          5      0   100%
# ptools/atom.py            159      0   100%
# ptools/attract.py          46     38    17%   23-33, 56-98
# ptools/forcefield.py      103     37    64%   140-193
# ptools/heligeom.py         26      0   100%
# ptools/io/__init__.py       2      0   100%
# ptools/io/attract.py       88      0   100%
# ptools/io/pdb.py           35      5    86%   56-57, 61-63
# ptools/pairlist.py         27      0   100%
# ptools/rigidbody.py        39      0   100%
# ptools/spatial.py         203      1    99%   401
# ptools/superpose.py       118     26    78%   54-55, 139-150, 153-164, 172, 190, 192, 196
# ptools/tables.py            1      0   100%
# -----------------------------------------------------
# TOTAL                     852    108    87%
