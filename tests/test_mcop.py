"""Tests for ptools.mcop"""

from tempfile import TemporaryFile
import unittest

from ptools.atom import AtomCollection, BaseAtom
from ptools.mcop import Mcop, McopRigid, McopRigidPDBError

from .testing.dummy import dummy_atomcollection, dummy_mcop, dummy_mcop_rigid
from .testing.io import mk_tmp_file, mk_pdb_models
from .testing.moreassert import assert_array_almost_equal, assert_array_not_almost_equal


class TestMcopPDBBuilder:
    """Helper class that builds dummy Mcop PDB files on demand.

    Model properties:

        - core : 11 atoms
        - regions : 2
            - region 1 : 4 atoms
            - region 2 : 3 atoms
    """

    @staticmethod
    def model_header(model_id: str) -> str:
        return f"MODEL     {model_id}"

    @staticmethod
    def model_footer() -> str:
        return "ENDMDL"

    @classmethod
    def _generic_region(cls, model_id: str, atoms: list[str]) -> str:
        return "\n".join([cls.model_header(model_id)] + atoms + [cls.model_footer()])

    @classmethod
    def core(cls):
        atoms = [
            "ATOM      1 CA   CYS     1      12.025  21.956  13.016    1   0.000 0 0",
            "ATOM      2 CSE  CYS     1      11.702  23.345  13.055    7   0.000 0 0",
            "ATOM      3 CA   GLY     2      12.408  20.728  16.555    1   0.000 0 0",
            "ATOM      4 CA   VAL     3      11.501  17.132  16.643    1   0.000 0 0",
            "ATOM      5 CSE  VAL     3      10.112  16.917  16.247   29   0.000 0 0",
            "ATOM      6 CA   PRO     4      14.215  14.528  17.062    1   0.000 0 0",
            "ATOM      7 CSE  PRO     4      14.229  15.158  18.286   22   0.000 0 0",
            "ATOM      8 CA   ALA     5      13.919  11.330  15.124    1   0.000 0 0",
            "ATOM      9 CSE  ALA     5      14.424  11.148  14.581    2   0.000 0 0",
            "ATOM     10 CA   ILE     6      15.540  10.094  18.303    1   0.000 0 0",
            "ATOM     11 CSE  ILE     6      17.227   9.373  18.051   14   0.000 0 0",
        ]
        return cls._generic_region(model_id="0", atoms=atoms)

    @classmethod
    def region1(cls):
        atoms = [
            "ATOM      1 CA   GLY   142      31.056  22.061  28.780    1   0.000 0 0",
            "ATOM      2 CA   LEU   143      32.876  20.969  32.008    1   0.000 0 0",
            "ATOM      3 CSE  LEU   143      34.750  21.987  32.256   15   0.000 0 0",
            "ATOM      4 CSE  LEU   143      34.750  21.987  32.256   15   0.000 0 0",
        ]
        return cls._generic_region(model_id="1  1", atoms=atoms)

    @classmethod
    def region2(cls):
        atoms = [
            "ATOM      1 CA   GLY   142      31.056  22.061  28.780    1   0.000 0 0",
            "ATOM      2 CA   LEU   143      32.876  20.969  32.008    1   0.000 0 0",
            "ATOM      3 CSE  LEU   143      34.750  21.987  32.256   15   0.000 0 0",
        ]
        return cls._generic_region(model_id="1  2", atoms=atoms)

    @classmethod
    def mcop_pdb(
        cls, has_core: bool = True, has_region1: bool = True, has_region2: bool = True
    ):
        regions = []
        if has_core:
            regions.append(cls.core())
        if has_region1:
            regions.append(cls.region1())
        if has_region2:
            regions.append(cls.region2())
        return "\n".join(regions)


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
        self._test_initialization_from_pdb(mcop)

    def test_from_pdb(self):
        with mk_pdb_models(n_models=3) as pdb:
            mcop = Mcop.from_pdb(pdb.name)
        self._test_initialization_from_pdb(mcop)

    def _test_initialization_from_pdb(self, mcop: Mcop):
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
        content = TestMcopPDBBuilder.mcop_pdb(has_core=False)
        with mk_tmp_file(content=content) as tmpfile:
            with self.assertRaisesRegex(
                McopRigidPDBError,
                "expecting core region first",
            ):
                McopRigid().read_pdb(tmpfile.name)

    def test_read_pdb_no_regions(self):
        content = TestMcopPDBBuilder.mcop_pdb(has_region1=False, has_region2=False)
        with mk_tmp_file(content=content) as tmpfile:
            with self.assertRaisesRegex(McopRigidPDBError, "no region found"):
                McopRigid().read_pdb(tmpfile.name)

    def test_read_pdb(self):
        content = TestMcopPDBBuilder.mcop_pdb()
        with mk_tmp_file(content=content) as tmpfile:
            rigid = McopRigid()
            rigid.read_pdb(tmpfile.name)
        self._test_initialization_from_pdb(rigid)

    def test_from_pdb(self):
        content = TestMcopPDBBuilder.mcop_pdb()
        with mk_tmp_file(content=content) as tmpfile:
            rigid = McopRigid.from_pdb(tmpfile.name)
        self._test_initialization_from_pdb(rigid)

    def _test_initialization_from_pdb(self, rigid: McopRigid):
        """Asserts initialization from PDB file went as expected."""
        self.assertEqual(len(rigid.core), 11)
        self.assertEqual(len(rigid.regions), 2)
        self.assertEqual(len(rigid.regions[0]), 4)
        self.assertEqual(len(rigid.regions[1]), 3)



def dummy_atomcollection(n_atoms: int = 10) -> AtomCollection:
    """Creates a dummy AtomCollection composed of `n_atoms` atoms."""
    return AtomCollection([BaseAtom(coords=(i, i, i)) for i in range(n_atoms)])
