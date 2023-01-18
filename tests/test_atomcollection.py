"""test_atomcollection - Tests for `ptools.atom.AtomCollection`."""

import unittest

import numpy as np

import ptools
from ptools import tables
from ptools.atomattrs import AtomAttrs
from ptools.atomcollection import Atom, AtomCollection

from .testing.moreassert import (
    assert_array_equal,
    assert_array_almost_equal,
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
        self.assertEqual(self.atoms.coordinates.shape, (self.n_atoms, 3))

    def test_initialization_empty(self):
        atoms = AtomCollection()
        self.assertEqual(len(atoms), 0)
        self.assertEqual(atoms.coordinates.shape, (0, 3))

    def _assert_copy_successful(self, thecopy: AtomCollection):
        # Check that both AtomCollections have the same dimensions.
        self.assertEqual(len(thecopy), self.n_atoms)
        self.assertEqual(thecopy.coordinates.shape, (self.n_atoms, 3))

        # Change parent AtomCollection coordinates and make sure it does not
        # affect the copy.
        ref_coords = self.atoms.coordinates.copy()
        self.atoms.coordinates.fill(0)
        assert_array_equal(thecopy.coordinates, ref_coords)

    def test_copy(self):
        thecopy = self.atoms.copy()
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
        assert_array_almost_equal(self.atoms.coordinates, ref_coords)

    def test_set_atom_coordinates_from_array(self):
        self.atoms.coordinates[0] = (42, 42, 42)
        assert_array_almost_equal(self.atoms[0].coordinates, (42, 42, 42))

    def test_set_atom_coordinates_from_atom(self):
        self.atoms[0].coordinates = (42, 42, 42)
        assert_array_almost_equal(self.atoms.coordinates[0], (42, 42, 42))

    def test_getitem_single_atom(self):
        atom = self.atoms[0]
        self.assertIsInstance(atom, Atom)
        assert_array_almost_equal(atom.coordinates, self.atoms.coordinates[0])

    def test_getitem_slice(self):
        self.assertIsInstance(self.atoms[:5], AtomCollection)

    def test_getitem_doesnt_make_copies(self):
        atom = self.atoms[0]
        atom.coordinates += 10
        assert_array_almost_equal(atom.coordinates, self.atoms.coordinates[0])

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

    def test_add(self):
        atoms2 = generate_dummy_atomcollection()
        n_final = len(self.atoms) + len(atoms2)
        all_atoms = self.atoms + atoms2
        self.assertEqual(len(all_atoms), n_final)
        self.assertEqual(all_atoms.coordinates.shape[0], n_final)

    def test_add_makes_copies(self):
        atoms2 = generate_dummy_atomcollection()
        assert_array_almost_equal(atoms2[-1].coordinates, [9, 9, 9])

        # Concatenates the two AtomCollections
        all_atoms = self.atoms + atoms2

        # Modifying the coordinates of the last atom in the second AtomCollection
        # should not affect the coordinates of the last atom in the concatenated
        # AtomCollection.
        atoms2[-1].coordinates = np.zeros(3)
        assert_array_almost_equal(atoms2[-1].coordinates, [0, 0, 0])
        assert_array_almost_equal(all_atoms[-1].coordinates, [9, 9, 9])

    def test_iadd(self):
        atoms = self.atoms.copy()
        atoms += self.atoms
        self.assertEqual(len(atoms), self.n_atoms * 2)
        self.assertEqual(atoms.coordinates.shape[0], self.n_atoms * 2)

    def test_masses(self):
        assert_array_equal(
            self.atoms.masses, np.full((self.n_atoms,), tables.masses["C"])
        )

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

    def test_set_chain(self):
        self.atoms.set_chain("A")
        self.assertEqual("".join(a.chain for a in self.atoms), "A" * len(self.atoms))
        self.atoms.set_chain("B")
        self.assertEqual("".join(a.chain for a in self.atoms), "B" * len(self.atoms))


class TestGuessAtomName(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(AtomAttrs.guess_element("CA"), "C")
        self.assertEqual(AtomAttrs.guess_element("NZ"), "N")

    def test_has_numeric(self):
        self.assertEqual(AtomAttrs.guess_element("CA1"), "C")

    def test_starts_with_numeric(self):
        self.assertEqual(AtomAttrs.guess_element("1CA"), "C")

    def test_has_only_numeric(self):
        self.assertEqual(AtomAttrs.guess_element("111"), "X")


class TestGuessAtomMass(unittest.TestCase):
    def test_basic(self):
        self.assertAlmostEqual(AtomAttrs.guess_mass("C"), 12.01100)
        self.assertAlmostEqual(AtomAttrs.guess_mass("H"), 1.00800)
        self.assertAlmostEqual(AtomAttrs.guess_mass("O"), 15.99900)
        self.assertAlmostEqual(AtomAttrs.guess_mass("N"), 14.00700)
        self.assertAlmostEqual(AtomAttrs.guess_mass("P"), 30.97400)
        self.assertAlmostEqual(AtomAttrs.guess_mass("S"), 32.06000)

    def test_element_does_not_exists(self):
        self.assertAlmostEqual(AtomAttrs.guess_mass("DUMMY"), 1.0)
