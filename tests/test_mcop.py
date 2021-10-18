"""Tests for ptools.mcop"""

import unittest

from ptools.atom import AtomCollection, BaseAtom
from ptools.mcop import Mcop, McopRigid, McopRigidPDBError

from .testing.dummy import dummy_atomcollection, dummy_mcop, dummy_mcop_rigid
from .testing.io import mk_tmp_file, mk_pdb_models
from .testing.moreassert import assert_array_almost_equal, assert_array_not_almost_equal



PDB_MCOPRIGID = {
    "core": """\
MODEL       0
ATOM      1 CA   CYS     1      12.025  21.956  13.016    1   0.000 0 0
ATOM      2 CSE  CYS     1      11.702  23.345  13.055    7   0.000 0 0
ATOM      3 CA   GLY     2      12.408  20.728  16.555    1   0.000 0 0
ATOM      4 CA   VAL     3      11.501  17.132  16.643    1   0.000 0 0
ATOM      5 CSE  VAL     3      10.112  16.917  16.247   29   0.000 0 0
ATOM      6 CA   PRO     4      14.215  14.528  17.062    1   0.000 0 0
ATOM      7 CSE  PRO     4      14.229  15.158  18.286   22   0.000 0 0
ATOM      8 CA   ALA     5      13.919  11.330  15.124    1   0.000 0 0
ATOM      9 CSE  ALA     5      14.424  11.148  14.581    2   0.000 0 0
ATOM     10 CA   ILE     6      15.540  10.094  18.303    1   0.000 0 0
ATOM     11 CSE  ILE     6      17.227   9.373  18.051   14   0.000 0 0
""",
    "region1": """\
MODEL       1  1
ATOM      1 CA   GLY   142      31.056  22.061  28.780    1   0.000 0 0
ATOM      2 CA   LEU   143      32.876  20.969  32.008    1   0.000 0 0
ATOM      3 CSE  LEU   143      34.750  21.987  32.256   15   0.000 0 0
ENDMDL    1  1
""",
    "region2": """\
MODEL       1  2
ATOM      1 CA   GLY   142      31.056  22.061  28.780    1   0.000 0 0
ATOM      2 CA   LEU   143      32.876  20.969  32.008    1   0.000 0 0
ATOM      3 CSE  LEU   143      34.750  21.987  32.256   15   0.000 0 0
ENDMDL    1  2
"""
}


class TestMcop(unittest.TestCase):
    """Tests for the Mcop class."""

    def test_add_copy(self):
        mcop = Mcop()
        self.assertEqual(len(mcop.copies), 0)
        mcop.add_copy(dummy_atomcollection())
        self.assertEqual(len(mcop.copies), 1)

    def test_add_copy_makes_copies(self):
        col = dummy_atomcollection()
        mcop = Mcop()
        mcop.add_copy(col)
        col.coords += 12  # after modification, mcop copy should be different.
        self.assertFalse(assert_array_almost_equal(col.coords, mcop.copies[0].coords))

    def test_size(self):
        mcop = Mcop()
        self.assertEqual(mcop.size(), 0)
        mcop.add_copy(dummy_atomcollection())
        self.assertEqual(mcop.size(), 1)

    def test_len(self):
        mcop = Mcop()
        self.assertEqual(len(mcop), 0)
        mcop.add_copy(dummy_atomcollection())
        self.assertEqual(len(mcop), 1)

    def test_getitem(self):
        mcop = Mcop()
        atoms = dummy_atomcollection()
        mcop.add_copy(atoms)
        atoms2 = mcop[0]
        assert_array_almost_equal(atoms.coords, atoms2.coords)

    def test_clear(self):
        mcop = Mcop()
        mcop.add_copy(dummy_atomcollection())
        self.assertEqual(len(mcop), 1)
        mcop.clear()
        self.assertEqual(len(mcop), 0)

    def test_read_pdb_mcop(self):
        mcop = Mcop()
        with mk_pdb_models(n_models=3) as pdb:
            mcop.read_pdb(pdb.name)
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
        self.rigid.add_region(dummy_atomcollection())
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

    def test_read_pdb_no_core_region_first(self):
        content = f"{PDB_MCOPRIGID['region1']}\n{PDB_MCOPRIGID['region2']}"
        with mk_tmp_file(content=content) as tmpfile:
            with self.assertRaisesRegex(McopRigidPDBError, "expecting core region first"):
                McopRigid().read_pdb(tmpfile.name)

    def test_read_pdb_no_regions(self):
        content = PDB_MCOPRIGID["core"]
        with mk_tmp_file(content=content) as tmpfile:
            with self.assertRaisesRegex(McopRigidPDBError, "no region found"):
                McopRigid().read_pdb(tmpfile.name)




    # def test_read_pdb(self):
    #     McopRigid().read_pdb(TEST_PDB_MCOPRIGID)
    #     self.assertEqual("not implemented", "coucou")




def dummy_atomcollection(n_atoms: int = 10) -> AtomCollection:
    """Creates a dummy AtomCollection composed of `n_atoms` atoms."""
    return AtomCollection([BaseAtom(coords=(i, i, i)) for i in range(n_atoms)])
