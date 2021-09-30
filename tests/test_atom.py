"""test_atom - Tests for `ptools.atom` module."""

import unittest

from ptools.atom import Atom, AtomCollection

from .testing.moreassert import (
    assert_array_almost_equal,
    assert_dummy_atom_initialization_ok,
    dummy_atom
)

class TestAtom(unittest.TestCase):
    def setUp(self):
        self.atoms = AtomCollection([dummy_atom()])

    def test_constructor_initializes_an_Atom_instance(self):
        self.assertIsInstance(self.atoms[0], Atom)

    def test_constructor(self):
        """Checks that when using copy constructor, all arguments are
        adequatly set and that atom coordinates are not a reference to
        the initial atom coordinates."""
        assert_dummy_atom_initialization_ok(self, self.atoms[0])

    def test_constructor_copies_coordinates(self):
        orig = dummy_atom()
        _ = AtomCollection([orig])
        atom = self.atoms[0]
        assert_array_almost_equal(atom.coords, (1, 2, 3))
        orig.coords = (0, 0, 0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))

    def test_mass_getter(self):
        self.assertEqual(self.atoms[0].mass, 12.011)

    # pylint guesses type wrong (self.atoms[0] is an Atom)
    # E1101: Instance of 'AtomCollection' has no 'mass' member
    # pylint: disable=E1101
    def test_mass_setter(self):
        self.atoms[0].mass = 42
        self.assertEqual(self.atoms[0].mass, 42)

    def test_equal(self):
        # Identical atoms from different AtomCollection instances should be evaluated equal.
        left = self.atoms[0]
        right = AtomCollection([left.copy()])[0]
        self.assertEqual(left, right)
