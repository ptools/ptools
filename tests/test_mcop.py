"""Tests for ptools.mcop"""

import unittest

from ptools.atom import AtomCollection, BaseAtom
from ptools.mcop import Mcop

from .testing.moreassert import assert_array_almost_equal
from . import TEST_PDB_3MODELS


class TestMcop(unittest.TestCase):
    """Tests for the Mcop class."""

    def test_add_copy(self):
        mcop = Mcop()
        self.assertEqual(len(mcop.copies), 0)
        mcop.add_copy(dummy_atom_collection())
        self.assertEqual(len(mcop.copies), 1)

    def test_add_copy_makes_copies(self):
        col = dummy_atom_collection()
        mcop = Mcop()
        mcop.add_copy(col)
        col.coords += 12  # after modification, mcop copy should be different.
        self.assertFalse(assert_array_almost_equal(col.coords, mcop.copies[0].coords))

    def test_size(self):
        mcop = Mcop()
        self.assertEqual(mcop.size(), 0)
        mcop.add_copy(dummy_atom_collection())
        self.assertEqual(mcop.size(), 1)

    def test_len(self):
        mcop = Mcop()
        self.assertEqual(len(mcop), 0)
        mcop.add_copy(dummy_atom_collection())
        self.assertEqual(len(mcop), 1)

    def test_getitem(self):
        mcop = Mcop()
        atoms = dummy_atom_collection()
        mcop.add_copy(atoms)
        atoms2 = mcop[0]
        assert_array_almost_equal(atoms.coords, atoms2.coords)

    def test_clear(self):
        mcop = Mcop()
        mcop.add_copy(dummy_atom_collection())
        self.assertEqual(len(mcop), 1)
        mcop.clear()
        self.assertEqual(len(mcop), 0)

    def test_read_pdb_mcop(self):
        mcop = Mcop()
        mcop.read_pdb(TEST_PDB_3MODELS)
        self.assertEqual(len(mcop), 3)




def dummy_atom_collection(n_atoms=10) -> AtomCollection:
    """Create a dummy AtomCollection composed of `n_atoms` atoms."""
    return AtomCollection([BaseAtom(coords=(i, i, i)) for i in range(n_atoms)])
