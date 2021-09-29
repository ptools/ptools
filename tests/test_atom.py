"""test_atom - Tests for `ptools.atom` module."""

import unittest

import ptools
from ptools.atom import BaseAtom, AtomCollection

from .testing.moreassert import assert_array_almost_equal


class TestAtom(unittest.TestCase):
    def setUp(self):
        orig = BaseAtom(
            name="CA",
            resname="ALA",
            chain="A",
            index=42,
            resid=17,
            charge=2.0,
            coords=(1, 2, 3),
        )
        self.collection = AtomCollection([orig])

    def test_constructor_initializes_an_Atom_instance(self):
        self.assertIsInstance(self.collection.atoms[0], ptools.atom.Atom)

    def test_constructor(self):
        """Checks that when using copy constructor, all arguments are
        adequatly set and that atom coordinates are not a reference to
        the initial atom coordinates."""
        atom = self.collection.atoms[0]
        self.assertEqual(atom.name, "CA")
        self.assertEqual(atom.resname, "ALA")
        self.assertEqual(atom.chain, "A")
        self.assertEqual(atom.index, 42)
        self.assertEqual(atom.resid, 17)
        self.assertEqual(atom.charge, 2.0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))

    def test_constructor_copies_coordinates(self):
        orig = BaseAtom(
            name="CA",
            resname="ALA",
            chain="A",
            index=42,
            resid=17,
            charge=2.0,
            coords=(1, 2, 3),
        )
        _ = AtomCollection([orig])
        atom = self.collection.atoms[0]
        assert_array_almost_equal(atom.coords, (1, 2, 3))
        orig.coords = (0, 0, 0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))

    def test_mass_getter(self):
        self.assertEqual(self.collection.atoms[0].mass, 12.011)

    def test_mass_setter(self):
        self.collection.atoms[0].mass = 42
        self.assertEqual(self.collection.atoms[0].mass, 42)

    def test_equal(self):
        # Identical atoms from different AtomCollection instances should be evaluated equal.
        left = self.collection[0]
        right = AtomCollection([left.copy()])[0]
        self.assertEqual(left, right)
