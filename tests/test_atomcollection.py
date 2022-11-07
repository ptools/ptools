"""test_atomcollection - Tests for `ptools.atom.AtomCollection`."""


import tempfile
import unittest

import numpy as np

import ptools
from ptools import tables
from ptools.atom import BaseAtom
from ptools.atomcollection import AtomCollection

from .testing.moreassert import (
    assert_array_equal,
    assert_array_almost_equal,
    assert_array_not_almost_equal,
)
from .testing.dummy import generate_dummy_atomcollection
from . import TEST_LIGAND


# Ignores R0904: Too many public methods
# pylint: disable=R0904
class TestAtomCollection(unittest.TestCase):
    def setUp(self):
        self.n_atoms: int = 10
        self.atoms: AtomCollection = generate_dummy_atomcollection(self.n_atoms)

    def test_initialization(self):
        self.assertEqual(len(self.atoms), self.n_atoms)
        self.assertEqual(self.atoms.coords.shape, (self.n_atoms, 3))

    def test_initialization_empty(self):
        atoms = AtomCollection()
        self.assertEqual(len(atoms), 0)
        self.assertEqual(atoms.coords.shape, (0, 3))

    def _assert_copy_successful(self, thecopy):
        # Check that both AtomCollections have the same dimensions.
        self.assertEqual(len(thecopy), self.n_atoms)
        self.assertEqual(thecopy.coords.shape, (self.n_atoms, 3))

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
        self.assertEqual(len(self.atoms), self.n_atoms)

    def test_size(self):
        self.assertEqual(self.atoms.size(), self.n_atoms)

    def test_isiterable(self):
        # If an exeception is raised here, AtomCollection is not iterable
        # which is not what we want.
        iter(self.atoms)

    def test_coordinates(self):
        ref_coords = np.array([[i, i, i] for i in range(self.n_atoms)], dtype=float)
        assert_array_almost_equal(self.atoms.coords, ref_coords)

    def test_set_atom_coordinates_from_array(self):
        self.atoms.coords[0] = (42, 42, 42)
        assert_array_almost_equal(self.atoms[0].coords, (42, 42, 42))

    def test_set_atom_coordinates_from_atom(self):
        self.atoms[0].coords = (42, 42, 42)
        assert_array_almost_equal(self.atoms.coords[0], (42, 42, 42))

    def test_getitem_single_atom(self):
        atom = self.atoms[0]
        self.assertIsInstance(atom, ptools.atom.Atom)
        assert_array_almost_equal(atom.coords, self.atoms.coords[0])

    def test_getitem_slice(self):
        self.assertIsInstance(self.atoms[:5], AtomCollection)

    def test_getitem_doesnt_make_copies(self):
        atom = self.atoms[0]
        atom.coords += 10
        assert_array_almost_equal(atom.coords, self.atoms.coords[0])

    def test_contains(self):
        atom = self.atoms[0]
        self.assertIn(atom, self.atoms)
        atom = self.atoms[0].copy()  # should return True as well
        self.assertIn(atom, self.atoms)

    def test_update_name(self):
        # Updating an atom's name should also update its type and mass.
        for i in range(self.n_atoms):
            self.atoms[i].name = "CA"
        self.assertEqual(
            "".join(at.element for at in self.atoms), "C" * len(self.atoms)
        )

        mass_ref = ptools.tables.masses["C"]
        assert_array_almost_equal(
            self.atoms.masses, np.ones(len(self.atoms)) * mass_ref
        )

    def test_centroid(self):
        center = self.atoms.centroid()
        assert_array_almost_equal(center, [4.5, 4.5, 4.5])

    def test_center_of_masses(self):
        # Masses are 1. COM should be same as center
        center = self.atoms.center_of_mass()
        assert_array_almost_equal(center, [4.5, 4.5, 4.5])

        # Changes atoms name, therefore element, therefore mass.
        for i in range(5):
            self.atoms[i].name = "CA"
        for i in range(5, 10):
            self.atoms[i].name = "NZ"
        self.atoms.guess_masses()

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
        center = self.atoms.centroid()
        self.atoms.translate(origin - center)
        assert_array_almost_equal(self.atoms.centroid(), (0, 0, 0))

    def test_translate_scalar(self):
        scalar = -4.5
        self.atoms.translate(scalar)
        assert_array_almost_equal(self.atoms.centroid(), (0, 0, 0))

    def test_center_without_weigths(self):
        origin = np.zeros(3)
        for origin in (np.zeros(3), np.ones(3)):
            assert_array_not_almost_equal(self.atoms.centroid(), origin)
            self.atoms.center_to_origin(origin, use_weights=False)
            assert_array_almost_equal(self.atoms.centroid(), origin)

    def test_center_with_weights(self):
        # Changes atom names, therefore elements, therefore masses.
        for i in range(5):
            self.atoms[i].name = "CA"
        for i in range(5, 10):
            self.atoms[i].name = "NZ"
        self.atoms.guess_masses()

        # This test should pass for a valid test of masses impact on AtomCollection.center()
        assert_array_not_almost_equal(
            self.atoms.centroid(), self.atoms.center_of_mass()
        )

        origin = np.zeros(3)
        for origin in (np.zeros(3), np.ones(3)):
            assert_array_not_almost_equal(self.atoms.center_of_mass(), origin)
            self.atoms.center_to_origin(origin, use_weights=True)
            assert_array_almost_equal(self.atoms.center_of_mass(), origin)

    def test_add(self):
        atoms2 = AtomCollection(
            [BaseAtom(coords=(i + 100, i, i)) for i in range(self.n_atoms)]
        )
        n_final = len(self.atoms) + len(atoms2)
        all_atoms = self.atoms + atoms2
        self.assertEqual(len(all_atoms), n_final)
        self.assertEqual(all_atoms.coords.shape[0], n_final)

    def test_add_makes_copies(self):
        atoms2 = AtomCollection(
            [BaseAtom(coords=(i + 100, i, i)) for i in range(self.n_atoms)]
        )
        assert_array_almost_equal(atoms2[-1].coords, [109, 9, 9])
        all_atoms = self.atoms + atoms2
        atoms2[-1].coords = np.zeros(3)

        assert_array_almost_equal(atoms2[-1].coords, [0, 0, 0])
        assert_array_almost_equal(all_atoms[-1].coords, [109, 9, 9])

    def test_iadd(self):
        atoms = self.atoms.copy()
        atoms += self.atoms
        self.assertEqual(len(atoms), self.n_atoms * 2)
        self.assertEqual(atoms.coords.shape[0], self.n_atoms * 2)

    def test_masses(self):
        assert_array_equal(
            self.atoms.masses, np.full((self.n_atoms,), tables.masses["C"])
        )

    @staticmethod
    def test_inertia_tensor():
        # Reference calculated with MDAnalysis 0.20.1:
        # >>> MDAnalysis.Universe("ligand.pdb").select_atoms("all").moment_of_inertia()
        atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)
        ref = [
            [3679339.47775172, 694837.16289436, -263651.10452372],
            [694837.16289436, 3803047.59374612, -194611.71739629],
            [-263651.10452372, -194611.71739629, 3425042.27240564],
        ]
        I = atoms.inertia_tensor()
        assert_array_almost_equal(I, ref, decimal=2)

    @staticmethod
    def test_principal_axes():
        atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)
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
        atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)
        sel = atoms.select_atom_type("CA")
        self.assertEqual(len(sel), 426)
        self.assertEqual([atom.name for atom in sel], ["CA"] * 426)

    def test_select_atom_types(self):
        atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)
        sel = atoms.select_atom_types(["CA", "CB"])
        self.assertEqual(len(sel), 426 + 64)
        self.assertEqual(
            [atom.name for atom in sel if atom.name.strip() == "CA"], ["CA"] * 426
        )
        self.assertEqual(
            [atom.name for atom in sel if atom.name.strip() == "CB"], ["CB"] * 64
        )

    def test_select_residue_range(self):
        atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)
        sel = atoms.select_residue_range(10, 20)
        self.assertEqual(len(sel), 23)
        self.assertTrue(all(10 <= atom.residue_index <= 20 for atom in sel))

    def test_select_chain(self):
        atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)
        sel = atoms.select_chain("B")
        self.assertEqual(len(sel), 974)
        self.assertTrue(all(atom.chain == "B" for atom in sel))

    def test_to_pdb(self):
        s = self.atoms.topdb()
        for i, line in enumerate(s.splitlines()):
            self.assertEqual(line, self.atoms[i].topdb())

    def test_writepdb(self):
        # Write PDB to temporary file.
        with tempfile.NamedTemporaryFile() as pdbfile:
            self.atoms.writepdb(pdbfile.name)

            # Read file back to check what's been written.
            with open(pdbfile.name, "rt", encoding="utf-8") as f:
                s = f.read().rstrip()

            self.assertEqual(s, self.atoms.topdb())

    def test_set_chain(self):
        self.atoms.set_chain("A")
        self.assertEqual("".join(a.chain for a in self.atoms), "A" * len(self.atoms))
        self.atoms.set_chain("B")
        self.assertEqual("".join(a.chain for a in self.atoms), "B" * len(self.atoms))


class TestGuessAtomName(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(ptools.atom.BaseAtom.guess_element("CA"), "C")
        self.assertEqual(ptools.atom.BaseAtom.guess_element("NZ"), "N")

    def test_has_numeric(self):
        self.assertEqual(ptools.atom.BaseAtom.guess_element("CA1"), "C")

    def test_starts_with_numeric(self):
        self.assertEqual(ptools.atom.BaseAtom.guess_element("1CA"), "C")

    def test_has_only_numeric(self):
        self.assertEqual(ptools.atom.BaseAtom.guess_element("111"), "X")


class TestGuessAtomMass(unittest.TestCase):
    def test_basic(self):
        self.assertAlmostEqual(ptools.atom.BaseAtom.guess_mass("C"), 12.01100)
        self.assertAlmostEqual(ptools.atom.BaseAtom.guess_mass("H"), 1.00800)
        self.assertAlmostEqual(ptools.atom.BaseAtom.guess_mass("O"), 15.99900)
        self.assertAlmostEqual(ptools.atom.BaseAtom.guess_mass("N"), 14.00700)
        self.assertAlmostEqual(ptools.atom.BaseAtom.guess_mass("P"), 30.97400)
        self.assertAlmostEqual(ptools.atom.BaseAtom.guess_mass("S"), 32.06000)

    def test_element_does_not_exists(self):
        self.assertAlmostEqual(ptools.atom.BaseAtom.guess_mass("DUMMY"), 1.0)
