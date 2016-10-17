
"""test_atom - Tests the `pyptools.atom` module."""

import unittest

from pyptools.atom import BaseAtom, AtomCollection

from . import assert_array_almost_equal


class TestAtom(unittest.TestCase):

    def test_empty_initializer(self):
        atom = BaseAtom()
        self.assertEqual(atom.name, 'XXX')
        self.assertEqual(atom.resname, 'XXX')
        self.assertEqual(atom.chain, 'X')
        self.assertEqual(atom.index, 0)
        self.assertEqual(atom.resid, 0)
        self.assertEqual(atom.charge, 0.0)
        assert_array_almost_equal(atom.coords, (0, 0, 0))

    def test_initialization(self):
        atom = BaseAtom(name='CA', resname='ALA', chain='A',
                        index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        self.assertEqual(atom.name, 'CA')
        self.assertEqual(atom.resname, 'ALA')
        self.assertEqual(atom.chain, 'A')
        self.assertEqual(atom.index, 42)
        self.assertEqual(atom.resid, 17)
        self.assertEqual(atom.charge, 2.0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))

    def test_initialize_with_bad_coordinates(self):
        err = 'atom coordinates but by a scalar or vector of shape 1 x 3 '
        with self.assertRaisesRegex(ValueError, err):
            BaseAtom(coords=(1, 2))

    def test_set_bad_coordinates(self):
        atom = BaseAtom()
        err = 'atom coordinates but by a scalar or vector of shape 1 x 3 '
        with self.assertRaisesRegex(ValueError, err):
            atom.coords = (1, 2)


class TestAtomCollection(unittest.TestCase):

    def test_initialization(self):
        atoms = [BaseAtom() for i in range(10)]
        atoms = AtomCollection(atoms)
        self.assertEqual(len(atoms), 10)
        # If an exeception is raised here, AtomCollection is not iterable
        # which is not what we want.
        iter(atoms)
