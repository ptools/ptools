"""test_atom - Tests for `ptools.atom.BaseAtom`."""

import unittest

from ptools.atom import BaseAtom
from ptools.array3d import Invalid3DArrayError
from .testing.moreassert import assert_array_almost_equal
from .testing.dummy import generate_dummy_atom


class TestBaseAtom(unittest.TestCase):
    def test_empty_initializer(self):
        atom = BaseAtom()
        self.assertEqual(atom.name, "XXX")
        self.assertEqual(atom.index, 0)
        self.assertEqual(atom.residue_name, "XXX")
        self.assertEqual(atom.residue_index, 0)
        self.assertEqual(atom.chain, "X")
        self.assertEqual(atom.charge, 0.0)
        assert_array_almost_equal(atom.coords, (0, 0, 0))

    def test_initialize_with_bad_coordinates(self):
        err = r"cannot initialize 3D-coordinates from array with shape \(2,\)"
        with self.assertRaisesRegex(Invalid3DArrayError, err):
            BaseAtom(coordinates=(1, 2))

    def test_set_bad_coordinates(self):
        atom = BaseAtom()
        err = r"cannot initialize 3D-coordinates from array with shape \(2,\)"
        with self.assertRaisesRegex(Invalid3DArrayError, err):
            atom.coords = (1, 2)

    # Ignores R0201: Method could be a function (no-self-use)
    # pylint: disable=R0201
    def test_set_coordinates(self):
        atom = BaseAtom()
        assert_array_almost_equal(atom.coords, (0, 0, 0))
        atom.coords = (1, 2, 3)
        assert_array_almost_equal(atom.coords, (1, 2, 3))

    def test_equals(self):
        self.assertEqual(BaseAtom(), BaseAtom())

    def test_equals_attributes_differ(self):
        left, right = BaseAtom(name="X"), BaseAtom(name="Y")
        self.assertNotEqual(left, right)

    def test_equals_coordinates_differ(self):
        left, right = BaseAtom(coordinates=(0, 0, 0)), BaseAtom(coordinates=(1, 1, 1))
        self.assertNotEqual(left, right)

    def test_copy(self):
        # Check that all arguments are adequatly set from original atom
        # and that atom coordinates are not a reference to
        # the initial atom coordinates.
        atom = generate_dummy_atom()
        original_coordinates = atom.coords.copy()
        atom_copy = atom.copy()
        self.assertNotEqual(
            atom_copy.coords.__array_interface__["data"][0],
            atom.coords.__array_interface__["data"][0],
        )
        assert_array_almost_equal(atom_copy.coords, original_coordinates)
        atom.coords = (0, 0, 0)
        assert_array_almost_equal(atom_copy.coords, original_coordinates)

    def test_name_setter(self):
        atom = BaseAtom()
        self.assertEqual(atom.element, "X")
        atom.name = "CA"
        self.assertEqual(atom.element, "C")

    def test_topdb(self):
        atom = generate_dummy_atom()
        reference_string = (
            "ATOM     42  CA  ALA A  17       "
            "1.000   2.000   3.000  1.00  0.00           "
            "C"
        )
        self.assertEqual(atom.topdb(), reference_string)

    def test_topdb_long_atomid(self):
        atom = generate_dummy_atom()
        atom.index = 110000
        reference_string = (
            "ATOM  1adb0  CA  ALA A  17       "
            "1.000   2.000   3.000  1.00  0.00           "
            "C"
        )
        self.assertEqual(atom.topdb(), reference_string)

    def test_topdb_long_resid(self):
        atom = generate_dummy_atom()
        atom.residue_index = 11000
        reference_string = (
            "ATOM     42  CA  ALA A2af8       "
            "1.000   2.000   3.000  1.00  0.00           "
            "C"
        )
        self.assertEqual(atom.topdb(), reference_string)

    def test_topdb_long_atom_name(self):
        atom = generate_dummy_atom()
        atom.name = "CA1"
        reference_string = (
            "ATOM     42  CA1 ALA A  17       "
            "1.000   2.000   3.000  1.00  0.00           "
            "C"
        )
        self.assertEqual(atom.topdb(), reference_string)
