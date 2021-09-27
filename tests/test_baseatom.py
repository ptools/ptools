"""test_atom - Tests for `ptools.atom.BaseAtom`."""

import unittest

from ptools.atom import BaseAtom
from .testing.moreassert import assert_array_almost_equal


class TestBaseAtom(unittest.TestCase):
    def test_empty_initializer(self):
        atom = BaseAtom()
        self.assertEqual(atom.name, "XXX")
        self.assertEqual(atom.resname, "XXX")
        self.assertEqual(atom.chain, "X")
        self.assertEqual(atom.index, 0)
        self.assertEqual(atom.resid, 0)
        self.assertEqual(atom.charge, 0.0)
        assert_array_almost_equal(atom.coords, (0, 0, 0))

    def test_initialization(self):
        atom = BaseAtom(
            name="CA",
            resname="ALA",
            chain="A",
            index=42,
            resid=17,
            charge=2.0,
            coords=(1, 2, 3),
        )
        self.assertEqual(atom.name, "CA")
        self.assertEqual(atom.resname, "ALA")
        self.assertEqual(atom.chain, "A")
        self.assertEqual(atom.index, 42)
        self.assertEqual(atom.resid, 17)
        self.assertEqual(atom.charge, 2.0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))

    def test_initialize_with_bad_coordinates(self):
        err = "3D coordinate array should be N x 3"
        with self.assertRaisesRegex(ValueError, err):
            BaseAtom(coords=(1, 2))

    def test_set_bad_coordinates(self):
        atom = BaseAtom()
        err = "3D coordinate array should be N x 3"
        with self.assertRaisesRegex(ValueError, err):
            atom.coords = (1, 2)

    def test_set_coordinates(self):
        atom = BaseAtom()
        assert_array_almost_equal(atom.coords, (0, 0, 0))
        atom.coords = (1, 2, 3)
        assert_array_almost_equal(atom.coords, (1, 2, 3))
    
    def test_repr(self):
        atom = BaseAtom()
        self.assertNotIn("_", repr(atom))

    def test_copy(self):
        # Check that all arguments are adequatly set from original atom
        # and that atom coordinates are not a reference to
        # the initial atom coordinates.
        atom = BaseAtom(
            name="CA",
            resname="ALA",
            chain="A",
            index=42,
            resid=17,
            charge=2.0,
            coords=(1, 2, 3),
        )
        atom_copy = atom.copy()
        self.assertEqual(atom_copy.name, "CA")
        self.assertEqual(atom_copy.resname, "ALA")
        self.assertEqual(atom_copy.chain, "A")
        self.assertEqual(atom_copy.index, 42)
        self.assertEqual(atom_copy.resid, 17)
        self.assertEqual(atom_copy.charge, 2.0)
        self.assertNotEqual(
            atom_copy.coords.__array_interface__["data"][0],
            atom.coords.__array_interface__["data"][0],
        )
        assert_array_almost_equal(atom_copy.coords, (1, 2, 3))
        atom.coords = (0, 0, 0)
        assert_array_almost_equal(atom_copy.coords, (1, 2, 3))

    def test_name_setter(self):
        atom = BaseAtom()
        self.assertEqual(atom.element, "X")
        atom.name = "CA"
        self.assertEqual(atom.element, "C")

    def test_topdb(self):
        atom = BaseAtom(
            name="CA",
            resname="ALA",
            chain="A",
            index=42,
            resid=17,
            charge=2.0,
            coords=(1, 2, 3),
        )
        reference_string = (
            "ATOM     42  CA  ALA A  17       "
            "1.000   2.000   3.000  1.00  0.00           "
            "C"
        )
        self.assertEqual(atom.topdb(), reference_string)

    def test_topdb_long_atomid(self):
        atom = BaseAtom(
            name="CA",
            resname="ALA",
            chain="A",
            index=110000,
            resid=17,
            charge=2.0,
            coords=(1, 2, 3),
        )
        reference_string = (
            "ATOM  1adb0  CA  ALA A  17       "
            "1.000   2.000   3.000  1.00  0.00           "
            "C"
        )
        self.assertEqual(atom.topdb(), reference_string)

    def test_topdb_long_resid(self):
        atom = BaseAtom(
            name="CA",
            resname="ALA",
            chain="A",
            index=42,
            resid=11000,
            charge=2.0,
            coords=(1, 2, 3),
        )
        reference_string = (
            "ATOM     42  CA  ALA A2af8       "
            "1.000   2.000   3.000  1.00  0.00           "
            "C"
        )
        self.assertEqual(atom.topdb(), reference_string)

    def test_topdb_long_atom_name(self):
        atom = BaseAtom(
            name="CA1",
            resname="ALA",
            chain="A",
            index=42,
            resid=17,
            charge=2.0,
            coords=(1, 2, 3),
        )
        reference_string = (
            "ATOM     42  CA1 ALA A  17       "
            "1.000   2.000   3.000  1.00  0.00           "
            "C"
        )
        self.assertEqual(atom.topdb(), reference_string)
