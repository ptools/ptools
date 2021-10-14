"""Tests for ptools.mcop"""

import unittest

from ptools.atom import AtomCollection, BaseAtom
from ptools.mcop import Mcop, McopRigid

from .testing.moreassert import assert_array_almost_equal, assert_array_not_almost_equal
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

    def test_attract_euler_rotate(self):
        """Tests that coordinates of each copy actually changed in the same fashion."""
        mcop = dummy_mcop()
        original_coordinates = mcop[0].coords.copy()
        mcop.attract_euler_rotate(10, 20, 30)
        current_coordinates = mcop[0].coords.copy()

        # Checks coordinates have changed.
        assert_array_not_almost_equal(current_coordinates, original_coordinates)

        # Checks all copies have the same coordinates.
        for kopy in mcop.copies[1:]:
            assert_array_almost_equal(kopy.coords, current_coordinates)



class TestMcopRigid(unittest.TestCase):

    def setUp(self):
        self.rigid = dummy_mcop_rigid()

    def test_add_region(self):
        n = len(self.rigid.regions)
        self.rigid.add_region(dummy_atom_collection())
        self.assertEqual(len(self.rigid.regions), n + 1)

    def test_attract_euler_rotate(self):
        """Tests that coordinates of the core and each copy of each region actually
        changed in the same fashion."""
        core_original = self.rigid.core.coords.copy()
        copy_original = self.rigid.regions[0][0].coords.copy()
        self.rigid.attract_euler_rotate(10, 20, 30)
        copy_rotated = self.rigid.regions[0][0].coords.copy()

        # Checks coordinates have changed.
        assert_array_not_almost_equal(self.rigid.core.coords, core_original)
        assert_array_not_almost_equal(self.rigid.regions[0][0].coords, copy_original)

        # Checks all copies have the same coordinates.
        for region in self.rigid.regions:
            for kopy in region.copies:
                assert_array_almost_equal(kopy.coords, copy_rotated)


def dummy_atom_collection(n_atoms: int = 10) -> AtomCollection:
    """Creates a dummy AtomCollection composed of `n_atoms` atoms."""
    return AtomCollection([BaseAtom(coords=(i, i, i)) for i in range(n_atoms)])


def dummy_mcop(n_copies: int = 10) -> Mcop:
    """Creates a dummpy Mcop composed of `n_copies` dummy AtomCollection instances."""
    mcop = Mcop()
    for _ in range(n_copies):
        mcop.add_copy(dummy_atom_collection())
    return mcop


def dummy_mcop_rigid(n_regions: int = 10) -> McopRigid:
    """Creates a dummpy McopRigid composed of `n_regions` that are dummy Mcop instances."""
    rigid = McopRigid()
    rigid.core = dummy_atom_collection()
    for _ in range(n_regions):
        rigid.add_region(dummy_mcop())
    return rigid