
"""test_atom - Tests for `ptools.atom` module."""

import unittest
import tempfile

import numpy as np

import ptools
from ptools.atom import (BaseAtom, Atom, AtomCollection, guess_atom_element,
                           guess_atom_mass)

from .testing.moreassert import assert_array_equal, assert_array_almost_equal

from . import TEST_LIGAND


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
        atom = BaseAtom(name="CA", resname="ALA", chain="A",
                        index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        self.assertEqual(atom.name, "CA")
        self.assertEqual(atom.resname, "ALA")
        self.assertEqual(atom.chain, "A")
        self.assertEqual(atom.index, 42)
        self.assertEqual(atom.resid, 17)
        self.assertEqual(atom.charge, 2.0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))

    def test_initialize_with_bad_coordinates(self):
        err = "3-d coordinates should be a scalar or 1 x 3 shaped-array"
        with self.assertRaisesRegex(ValueError, err):
            BaseAtom(coords=(1, 2))

    def test_set_bad_coordinates(self):
        atom = BaseAtom()
        err = "3-d coordinates should be a scalar or 1 x 3 shaped-array"
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
        atom = BaseAtom(name="CA", resname="ALA", chain="A",
                        index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        atom_copy = atom.copy()
        self.assertEqual(atom_copy.name, "CA")
        self.assertEqual(atom_copy.resname, "ALA")
        self.assertEqual(atom_copy.chain, "A")
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
        parent = BaseAtom(name="CA", resname="ALA", chain="A",
                          index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        atom = BaseAtom(orig=parent)
        self.assertEqual(atom.name, "CA")
        self.assertEqual(atom.resname, "ALA")
        self.assertEqual(atom.chain, "A")
        self.assertEqual(atom.index, 42)
        self.assertEqual(atom.resid, 17)
        self.assertEqual(atom.charge, 2.0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))
        parent.coords = (0, 0, 0)
        assert_array_almost_equal(atom.coords, (1, 2, 3))

    def test_topdb(self):
        atom = BaseAtom(name="CA", resname="ALA", chain="A",
                        index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        reference_string = ("ATOM     42  CA  ALA A  17       "
                            "1.000   2.000   3.000  1.00  0.00           "
                            "C")
        self.assertEqual(atom.topdb(), reference_string)

    def test_topdb_long_atomid(self):
        atom = BaseAtom(name="CA", resname="ALA", chain="A",
                        index=110000, resid=17, charge=2.0, coords=(1, 2, 3))
        reference_string = ("ATOM  1adb0  CA  ALA A  17       "
                            "1.000   2.000   3.000  1.00  0.00           "
                            "C")
        self.assertEqual(atom.topdb(), reference_string)

    def test_topdb_long_resid(self):
        atom = BaseAtom(name="CA", resname="ALA", chain="A",
                        index=42, resid=11000, charge=2.0, coords=(1, 2, 3))
        reference_string = ("ATOM     42  CA  ALA A2af8       "
                            "1.000   2.000   3.000  1.00  0.00           "
                            "C")
        self.assertEqual(atom.topdb(), reference_string)

    def test_topdb_long_atom_name(self):
        atom = BaseAtom(name="CA1", resname="ALA", chain="A",
                        index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        reference_string = ("ATOM     42  CA1 ALA A  17       "
                            "1.000   2.000   3.000  1.00  0.00           "
                            "C")
        self.assertEqual(atom.topdb(), reference_string)


class TestAtom(unittest.TestCase):

    def test_constructor(self):
        # Check that when using copy constructor, all arguments are
        # adequatly set and that atom coordinates are not a reference to
        # the initial atom coordinates.
        orig = BaseAtom(name="CA", resname="ALA", chain="A",
                        index=42, resid=17, charge=2.0, coords=(1, 2, 3))
        collection = AtomCollection([orig])
        atom = collection.atoms[0]
        self.assertEqual(atom.name, "CA")
        self.assertEqual(atom.resname, "ALA")
        self.assertEqual(atom.chain, "A")
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

    def test_initialization(self):
        self.assertEqual(len(self.atoms.atoms), 10)
        self.assertEqual(self.atoms.coords.shape, (10, 3))

    def test_initialization_empty(self):
        atoms = AtomCollection()
        self.assertEqual(len(atoms.atoms), 0)
        self.assertEqual(atoms.coords.shape, (0, 3))

    def _assert_copy_successful(self, thecopy):
        # Check that both AtomCollections have the same dimensions.
        self.assertEqual(len(thecopy), len(self.atoms))
        self.assertEqual(thecopy.coords.shape, (10, 3))

        # Change parent AtomCollection coordinates and make sure it does not
        # affect the copy.
        ref_coords = self.atoms.coords.copy()
        self.atoms.coords.fill(0)
        assert_array_equal(thecopy.coords, ref_coords)

    def test_copy_initialization1(self):
        thecopy = self.atoms.copy()
        self._assert_copy_successful(thecopy)

    def test_copy_initialization2(self):
        thecopy = AtomCollection.copy(self.atoms)
        self._assert_copy_successful(thecopy)

    def test_len(self):
        self.assertEqual(len(self.atoms), 10)

    def test_size(self):
        self.assertEqual(self.atoms.size(), 10)

    def test_isiterable(self):
        # If an exeception is raised here, AtomCollection is not iterable
        # which is not what we want.
        iter(self.atoms)

    def test_coordinates(self):
        ref_coords = np.array([[i, i, i] for i in range(10)], dtype=float)
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

    def test_update_name(self):
        # Updating an atom's name should also update its type and mass.
        for i in range(10):
            self.atoms[i].name = "CA"
        self.assertEqual("".join(at.element for at in self.atoms),
                         "C" * len(self.atoms))

        mass_ref = ptools.tables.masses["C"]
        assert_array_almost_equal(self.atoms.masses(),
                                  np.ones(len(self.atoms)) * mass_ref)

    def test_center(self):
        center = self.atoms.center()
        assert_array_almost_equal(center, [4.5, 4.5, 4.5])

    def test_center_of_masses(self):
        # Masses are 1. COM should be same as center
        center = self.atoms.center_of_mass()
        assert_array_almost_equal(center, [4.5, 4.5, 4.5])

        # Set masses to 12 (in principle, updating an atom's name
        # should also update its type).
        for i in range(5):
            self.atoms[i].name = "CA"
        for i in range(5, 10):
            self.atoms[i].name = "NZ"
        center = self.atoms.center_of_mass()
        assert_array_almost_equal(center, [4.69179, 4.69179, 4.69179])

    def test_radius_of_gyration(self):
        rgyr = self.atoms.radius_of_gyration()
        # Reference value calculated with VMD.
        self.assertAlmostEqual(rgyr, 4.9749369621276855, places=6)

    def test_translate(self):
        # Translate is a method herited from `spatial.SpatialObject`.
        # Basically is should work on any child class.
        origin = (0, 0, 0)
        center = self.atoms.center()
        self.atoms.translate(origin - center)
        assert_array_almost_equal(self.atoms.center(), (0, 0, 0))

    def test_translate_scalar(self):
        scalar = -4.5
        self.atoms.translate(scalar)
        assert_array_almost_equal(self.atoms.center(), (0, 0, 0))

    def test_add(self):
        atoms2 = AtomCollection([BaseAtom(coords=(i+100, i, i))
                                 for i in range(10)])
        all_atoms = self.atoms + atoms2
        N = len(self.atoms) + len(atoms2)
        self.assertEqual(len(all_atoms), N)

    def test_masses(self):
        # atoms name "XXX" should weight 0
        assert_array_equal(self.atoms.masses(), np.ones(10))

    def test_tensor_of_inertia(self):
        # Calculated with MDAnalysis 0.20.1:
        # >>> MDAnalysis.Universe("ligand.pdb").select_atoms("all").moment_of_inertia()
        atoms = ptools.io.read_pdb(TEST_LIGAND)
        ref = [[3679339.47775172,  694837.16289436, -263651.10452372],
               [ 694837.16289436, 3803047.59374612, -194611.71739629],
               [-263651.10452372, -194611.71739629, 3425042.27240564]]
        I = atoms.tensor_of_inertia()
        assert_array_almost_equal(I, ref, decimal=2)

    def test_principal_axes(self):
        atoms = ptools.io.read_pdb(TEST_LIGAND)
        # Calculated with MDAnalysis 0.20.1:
        # >>> MDAnalysis.Universe("ligand.pdb").select_atoms("all").principal_axes()
        ref = [[ 0.65682984,  0.70033642, -0.27946997],
               [ 0.04052252,  0.33731064,  0.94052084],
               [ 0.7529492 , -0.62908698,  0.1931763 ]]
        I = atoms.principal_axes()
        assert_array_almost_equal(I, ref)

    def test_to_pdb(self):
        s = self.atoms.topdb()
        for i, line in enumerate(s.splitlines()):
            self.assertEqual(line, self.atoms[i].topdb())

    def test_writepdb(self):
        # Write PDB to temporary file.
        f = tempfile.NamedTemporaryFile()
        self.atoms.writepdb(f.name)

        # Read file back to check what's been written.
        with open(f.name, "rt") as f:
            s = f.read().rstrip()

        self.assertEqual(s, self.atoms.topdb())

        # Close (and delete) temporary file.
        f.close()

    def test_set_chain(self):
        self.atoms.set_chain("A")
        self.assertEqual("".join(a.chain for a in self.atoms), "A" * len(self.atoms))
        self.atoms.set_chain("B")
        self.assertEqual("".join(a.chain for a in self.atoms), "B" * len(self.atoms))


class TestGuessAtomName(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(guess_atom_element("CA"), "C")
        self.assertEqual(guess_atom_element("NZ"), "N")

    def test_has_numeric(self):
        self.assertEqual(guess_atom_element("CA1"), "C")

    def test_starts_with_numeric(self):
        self.assertEqual(guess_atom_element("1CA"), "C")

    def test_has_only_numeric(self):
        self.assertEqual(guess_atom_element("111"), "X")


class TestGuessAtomMass(unittest.TestCase):
    def test_basic(self):
        self.assertAlmostEqual(guess_atom_mass("C"), 12.01100)
        self.assertAlmostEqual(guess_atom_mass("H"), 1.00800)
        self.assertAlmostEqual(guess_atom_mass("O"), 15.99900)
        self.assertAlmostEqual(guess_atom_mass("N"), 14.00700)
        self.assertAlmostEqual(guess_atom_mass("P"), 30.97400)
        self.assertAlmostEqual(guess_atom_mass("S"), 32.06000)

    def test_element_does_not_exists(self):
        self.assertAlmostEqual(guess_atom_mass("DUMMY"), 1.0)
