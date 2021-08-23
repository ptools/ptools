"""test_atomcollection - Tests for `ptools.atom.AtomCollection`."""

import unittest
import tempfile

import numpy as np

import ptools
from ptools.atom import AtomCollection, BaseAtom

from .testing.moreassert import assert_array_equal, assert_array_almost_equal
from . import TEST_LIGAND


class TestAtomCollection(unittest.TestCase):
    def setUp(self):
        # Create an AtomCollection with 10 atoms.
        # Atom coordinates are [(0, 0, 0), (1, 1, 1), ..., (9, 9, 9)].
        self.atoms = AtomCollection([BaseAtom(coords=(i, i, i)) for i in range(10)])

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

    def test_repr(self):
        s = repr(self.atoms)
        self.assertIn("AtomCollection", s)
        self.assertIn(str(len(self.atoms)), s)

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
        self.assertIsInstance(atom, ptools.atom.Atom)
        assert_array_almost_equal(atom.coords, self.atoms.coords[0])
        atom.coords = [42, 17, 323]
        assert_array_almost_equal(atom.coords, self.atoms.coords[0])

    def test_update_name(self):
        # Updating an atom's name should also update its type and mass.
        for i in range(10):
            self.atoms[i].name = "CA"
        self.assertEqual(
            "".join(at.element for at in self.atoms), "C" * len(self.atoms)
        )

        mass_ref = ptools.tables.masses["C"]
        assert_array_almost_equal(
            self.atoms.masses, np.ones(len(self.atoms)) * mass_ref
        )

    def test_center(self):
        center = self.atoms.center()
        assert_array_almost_equal(center, [4.5, 4.5, 4.5])

    def test_center_of_masses(self):
        # Masses are 1. COM should be same as center
        center = self.atoms.center_of_mass()
        assert_array_almost_equal(center, [4.5, 4.5, 4.5])

        # Update atoms name, therefore element, therefore mass.
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
        atoms2 = AtomCollection([BaseAtom(coords=(i + 100, i, i)) for i in range(10)])
        all_atoms = self.atoms + atoms2
        N = len(self.atoms) + len(atoms2)
        self.assertEqual(len(all_atoms), N)

    def test_add_makes_copies(self):
        atoms2 = AtomCollection([BaseAtom(coords=(i + 100, i, i)) for i in range(10)])
        assert_array_almost_equal(atoms2[-1].coords, [109, 9, 9])
        all_atoms = self.atoms + atoms2
        atoms2[-1].coords = np.zeros(3)

        assert_array_almost_equal(atoms2[-1].coords, [0, 0, 0])
        assert_array_almost_equal(all_atoms[-1].coords, [109, 9, 9])

    def test_masses(self):
        # atoms name "XXX" should weight 0
        assert_array_equal(self.atoms.masses, np.ones(10))

    def test_tensor_of_inertia_accurate(self):
        # Reference calculated with MDAnalysis 0.20.1:
        # >>> MDAnalysis.Universe("ligand.pdb").select_atoms("all").moment_of_inertia()
        atoms = ptools.io.read_pdb(TEST_LIGAND)
        ref = [
            [3679339.47775172, 694837.16289436, -263651.10452372],
            [694837.16289436, 3803047.59374612, -194611.71739629],
            [-263651.10452372, -194611.71739629, 3425042.27240564],
        ]
        I = atoms.tensor_of_inertia(method="accurate")
        assert_array_almost_equal(I, ref, decimal=2)

    def test_tensor_of_inertia_fast(self):
        # Reference calculated with ptools-python d0b41dc (this is actually a
        # non-regression test).
        atoms = ptools.io.read_pdb(TEST_LIGAND)
        ref = [
            [147200.90378622, -57794.00682372, 21649.47873511],
            [-57794.00682372, 136694.86723672, 16076.77968389],
            [21649.47873511, 16076.77968389, 168066.83235232],
        ]
        I = atoms.tensor_of_inertia(method="fast")
        assert_array_almost_equal(I, ref, decimal=8)

    def test_principal_axes(self):
        atoms = ptools.io.read_pdb(TEST_LIGAND)
        # Calculated with MDAnalysis 0.20.1:
        # >>> MDAnalysis.Universe("ligand.pdb").select_atoms("all").principal_axes()
        ref = [
            [0.65682984, 0.70033642, -0.27946997],
            [0.04052252, 0.33731064, 0.94052084],
            [0.7529492, -0.62908698, 0.1931763],
        ]
        I = atoms.principal_axes()
        assert_array_almost_equal(I, ref)

    def test_select_atom_type(self):
        atoms = ptools.io.read_pdb(TEST_LIGAND)
        sel = atoms.select_atom_type("CA")
        self.assertEqual(len(sel), 426)
        self.assertEqual([atom.name for atom in sel], ["CA"] * 426)

    def test_select_residue_range(self):
        atoms = ptools.io.read_pdb(TEST_LIGAND)
        sel = atoms.select_residue_range(10, 20)
        self.assertEqual(len(sel), 23)
        self.assertTrue(all(10 <= atom.resid <= 20 for atom in sel))

    def test_select_chain(self):
        atoms = ptools.io.read_pdb(TEST_LIGAND)
        sel = atoms.select_chain("B")
        self.assertEqual(len(sel), 974)
        self.assertTrue(all(atom.chain == "B" for atom in sel))

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
        self.assertEqual(ptools.atom.guess_atom_element("CA"), "C")
        self.assertEqual(ptools.atom.guess_atom_element("NZ"), "N")

    def test_has_numeric(self):
        self.assertEqual(ptools.atom.guess_atom_element("CA1"), "C")

    def test_starts_with_numeric(self):
        self.assertEqual(ptools.atom.guess_atom_element("1CA"), "C")

    def test_has_only_numeric(self):
        self.assertEqual(ptools.atom.guess_atom_element("111"), "X")


class TestGuessAtomMass(unittest.TestCase):
    def test_basic(self):
        self.assertAlmostEqual(ptools.atom.guess_atom_mass("C"), 12.01100)
        self.assertAlmostEqual(ptools.atom.guess_atom_mass("H"), 1.00800)
        self.assertAlmostEqual(ptools.atom.guess_atom_mass("O"), 15.99900)
        self.assertAlmostEqual(ptools.atom.guess_atom_mass("N"), 14.00700)
        self.assertAlmostEqual(ptools.atom.guess_atom_mass("P"), 30.97400)
        self.assertAlmostEqual(ptools.atom.guess_atom_mass("S"), 32.06000)

    def test_element_does_not_exists(self):
        self.assertAlmostEqual(ptools.atom.guess_atom_mass("DUMMY"), 1.0)
