
"""test_atom - Tests for `pyptools.atom` module."""

import unittest

import numpy

from pyptools.atom import BaseAtom, Atom, AtomCollection

from .testing.moreassert import assert_array_equal, assert_array_almost_equal


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
        # Check that when using copy constructor, all arguments are
        # adequatly set and that atom coordinates are not a reference to
        # the initial atom coordinates.
        orig = BaseAtom(name='CA', resname='ALA', chain='A',
                        index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        collection = AtomCollection([orig])
        atom = collection.atoms[0]
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
        # Create an AtomCollection with 10 atoms.
        # Atom coordinates are [(0, 0, 0), (1, 1, 1), ..., (9, 9, 9)].
        self.atoms = AtomCollection([BaseAtom(coords=(i, i, i))
                                     for i in range(10)])

    def test_constructor(self):
        self.assertEqual(len(self.atoms.atoms), 10)
        self.assertEqual(self.atoms.coords.shape, (10, 3))

    def _assert_copy_successful(self, thecopy):
        # Check that both AtomCollections have the same dimensions.
        self.assertEqual(len(thecopy), len(self.atoms))
        self.assertEqual(thecopy.coords.shape, (10, 3))

        # Change parent AtomCollection coordinates and make sure it does not
        # affect the copy.
        ref_coords = self.atoms.coords.copy()
        self.atoms.coords.fill(0)
        assert_array_equal(thecopy.coords, ref_coords)

    def test_copy_constructor1(self):
        thecopy = self.atoms.copy()
        self._assert_copy_successful(thecopy)

    def test_copy_constructor2(self):
        thecopy = AtomCollection.copy(self.atoms)
        self._assert_copy_successful(thecopy)

    def test_len(self):
        self.assertEqual(len(self.atoms), 10)

    def test_size(self):
        self.assertEqual(len(self.atoms), 10)

    def test_isiterable(self):
        # If an exeception is raised here, AtomCollection is not iterable
        # which is not what we want.
        iter(self.atoms)

    def test_coordinates(self):
        ref_coords = numpy.array([[i, i, i] for i in range(10)], dtype=float)
        assert_array_almost_equal(self.atoms.coords, ref_coords)

    def test_set_atom_coordinates_from_array(self):
        self.atoms.coords[0] = (42, 42, 42)
        assert_array_almost_equal(self.atoms.atoms[0].coords, (42, 42, 42))

    def test_set_atom_coordinates_from_atom(self):
        self.atoms.atoms[0].coords = (42, 42, 42)
        assert_array_almost_equal(self.atoms.coords[0], (42, 42, 42))

    def test_getitem(self):
        atom = self.atoms[0]
        self.assertTrue(isinstance(atom, Atom))
        assert_array_almost_equal(atom.coords, self.atoms.coords[0])
        atom.coords = [42, 17, 323]
        assert_array_almost_equal(atom.coords, self.atoms.coords[0])

    def test_get_center(self):
        center = self.atoms.get_center()
        assert_array_almost_equal(center, [4.5, 4.5, 4.5])

    def test_get_radius_of_gyration(self):
        rgyr = self.atoms.get_radius_of_gyration()
        # Reference value calculated with VMD.
        self.assertAlmostEqual(rgyr, 4.9749369621276855, places=6)

    def test_translate(self):
        # Translate is a method herited from `spatial.SpatialObject`.
        # Basically is should work on any child class.
        origin = (0, 0, 0)
        center = self.atoms.get_center()
        self.atoms.translate(origin - center)
        assert_array_almost_equal(self.atoms.get_center(), (0, 0, 0))

    def test_translate_scalar(self):
        scalar = -4.5
        self.atoms.translate(scalar)
        assert_array_almost_equal(self.atoms.get_center(), (0, 0, 0))
