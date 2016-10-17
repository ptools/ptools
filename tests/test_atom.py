
"""test_atom - Tests the `pyptools.atom` module."""

import unittest

from pyptools.atom import BaseAtom, Atom, AtomCollection

from . import assert_array_almost_equal


class TestBaseAtom(unittest.TestCase):

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
        err = '3-d coordinates should be a scalar or 1 x 3 shaped-array'
        with self.assertRaisesRegex(ValueError, err):
            BaseAtom(coords=(1, 2))

    def test_set_bad_coordinates(self):
        atom = BaseAtom()
        err = '3-d coordinates should be a scalar or 1 x 3 shaped-array'
        with self.assertRaisesRegex(ValueError, err):
            atom.coords = (1, 2)

    def test_set_coordinates(self):
        atom = BaseAtom()
        assert_array_almost_equal(atom.coords, (0, 0, 0))
        atom.coords = (1, 2, 3)
        assert_array_almost_equal(atom.coords, (1, 2, 3))

    def test_copy(self):
        # Check that all arguments are adequatly set from original atom
        # and that atom coordinates are not a reference to
        # the initial atom coordinates.
        atom = BaseAtom(name='CA', resname='ALA', chain='A',
                        index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        atom_copy = atom.copy()
        self.assertEqual(atom_copy.name, 'CA')
        self.assertEqual(atom_copy.resname, 'ALA')
        self.assertEqual(atom_copy.chain, 'A')
        self.assertEqual(atom_copy.index, 42)
        self.assertEqual(atom_copy.resid, 17)
        self.assertEqual(atom_copy.charge, 2.0)
        assert_array_almost_equal(atom_copy.coords, (1, 2, 3))
        atom.coords = (0, 0, 0)
        assert_array_almost_equal(atom_copy.coords, (1, 2, 3))

    def test_copy_constructor(self):
        # Check that when using copy constructor, all arguments are
        # adequatly set and that atom coordinates are not a reference to
        # the initial atom coordinates.
        parent = BaseAtom(name='CA', resname='ALA', chain='A',
                          index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        atom = BaseAtom(orig=parent)
        self.assertEqual(atom.name, 'CA')
        self.assertEqual(atom.resname, 'ALA')
        self.assertEqual(atom.chain, 'A')
        self.assertEqual(atom.index, 42)
        self.assertEqual(atom.resid, 17)
        self.assertEqual(atom.charge, 2.0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))
        parent.coords = (0, 0, 0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))


class TestAtom(unittest.TestCase):

    def test_constructor(self):
        # Check that all arguments are adequatly set from original atom
        # and that atom coordinates are not a reference to
        # the initial atom coordinates.
        orig = BaseAtom(name='CA', resname='ALA', chain='A',
                        index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        atom = Atom(orig)
        self.assertEqual(atom.name, 'CA')
        self.assertEqual(atom.resname, 'ALA')
        self.assertEqual(atom.chain, 'A')
        self.assertEqual(atom.index, 42)
        self.assertEqual(atom.resid, 17)
        self.assertEqual(atom.charge, 2.0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))
        orig.coords = (0, 0, 0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))


class TestAtomCollection(unittest.TestCase):

    def setUp(self):
        self.atoms = AtomCollection([BaseAtom() for i in range(10)])

    def test_initialization(self):
        self.assertEqual(len(self.atoms.atoms), 10)

    def test_len(self):
        self.assertEqual(len(self.atoms), 10)

    def test_isiterable(self):
        # If an exeception is raised here, AtomCollection is not iterable
        # which is not what we want.
        iter(self.atoms)
