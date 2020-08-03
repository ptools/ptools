
"""test_io - Tests for `ptools.io` module."""

import unittest

import ptools.io as io

from . import TEST_PDB, TEST_PDB_ATOM_NAMES
from .testing.io import random_filename


class TestIO(unittest.TestCase):

    def test_read_pdb_does_not_exist(self):
        filename = random_filename()
        with self.assertRaises(FileNotFoundError):
            io.read_pdb(filename)

    def test_read_pdb(self):
        atoms = io.read_pdb(TEST_PDB)
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
