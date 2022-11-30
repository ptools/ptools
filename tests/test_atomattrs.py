"""test_atomattrs - Tests for `ptools.atom.AtomAttrs`."""

import unittest

from ptools.atomattrs import AtomAttrs
from ptools.array3d import Invalid3DArrayError
from .testing.moreassert import assert_array_almost_equal, assert_array_not_almost_equal
from .testing.dummy import generate_dummy_atom


class TestAtomProperties(unittest.TestCase):
    def test_empty_initializer(self):
        atom = AtomAttrs()
        self.assertEqual(atom.name, "XXX")
        self.assertEqual(atom.index, 0)
        self.assertEqual(atom.residue_name, "XXX")
        self.assertEqual(atom.residue_index, 0)
        self.assertEqual(atom.chain, "X")
        self.assertEqual(atom.charge, 0.0)
        assert_array_almost_equal(atom.coordinates, (0, 0, 0))

    def test_initialize_with_bad_coordinates(self):
        err = r"cannot initialize 3D-coordinates from array with shape \(2,\)"
        with self.assertRaisesRegex(Invalid3DArrayError, err):
            AtomAttrs(coordinates=(1, 2))

    def test_set_bad_coordinates(self):
        atom = AtomAttrs()
        err = r"cannot initialize 3D-coordinates from array with shape \(2,\)"
        with self.assertRaisesRegex(Invalid3DArrayError, err):
            atom.coordinates = (1, 2)

    def test_set_coordinates(self):
        atom = AtomAttrs()
        assert_array_almost_equal(atom.coordinates, (0, 0, 0))
        atom.coordinates = (1, 2, 3)
        assert_array_almost_equal(atom.coordinates, (1, 2, 3))

    def test_equals(self):
        self.assertEqual(AtomAttrs(), AtomAttrs())

    def test_equals_attributes_differ(self):
        left, right = AtomAttrs(name="X"), AtomAttrs(name="Y")
        self.assertNotEqual(left, right)

    def test_equals_coordinates_differ(self):
        left, right = AtomAttrs(coordinates=(0, 0, 0)), AtomAttrs(coordinates=(1, 1, 1))
        self.assertNotEqual(left, right)

    def test_copy(self):
        atom = generate_dummy_atom()
        copy = atom.copy()

        # Checks equality in the standard fashion for all attributes expect ``coordinates``.
        for attr in atom.__dict__.keys() - ("coordinates",):
            self.assertEqual(getattr(atom, attr), getattr(copy, attr))

        # Checks numpy array of coordinates is actually a copy.
        self.assertNotEqual(id(atom.coordinates), id(copy.coordinates))
        assert_array_almost_equal(atom.coordinates, copy.coordinates)

        source_coordinates = atom.coordinates.copy()
        atom.coordinates += 12

        # copy coordinates should stay unchanged.
        assert_array_almost_equal(copy.coordinates, source_coordinates)

        # ...while original atom coordinates should have changed.
        assert_array_not_almost_equal(atom.coordinates, source_coordinates)

    def test_name_setter(self):
        atom = AtomAttrs()
        self.assertEqual(atom.element, "X")
        atom.name = "CA"
        self.assertEqual(atom.element, "C")
