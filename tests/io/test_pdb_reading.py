"""Tests for PDB reading facilities."""

import unittest

import ptools.io.pdb as pdb
from ptools.atom import AtomCollection

from ..testing.io import (
    mk_pdb_models,
    mk_pdb_no_model,
    mk_pdb_10_atoms,
    random_filename,
    TestPDBBuilder,
)


class TestPDBIO(unittest.TestCase):
    def test_read_pdb_does_not_exist(self):
        filename = random_filename()
        with self.assertRaises(FileNotFoundError):
            pdb.read_pdb(filename)

    def test_is_atom_line(self):
        self.assertTrue(pdb.PDBLine("ATOM  <rest_of_the_line>").is_atom())
        self.assertTrue(pdb.PDBLine("HETATM<rest_of_the_line>").is_atom())
        self.assertFalse(pdb.PDBLine("ATOMI").is_atom())
        self.assertFalse(pdb.PDBLine("ANISOU").is_atom())

    def test_read_pdb(self):
        with mk_pdb_no_model() as pdb_file:
            atoms = pdb.read_pdb(pdb_file.name)
        self.assertEqual(len(atoms), 10)
        for i in range(10):
            atom = atoms[i]
            self.assertEqual(atom.index, i + 1)
            self.assertEqual(atom.name, TestPDBBuilder.atom_names()[i])
            self.assertEqual(atom.resid, 1)
            self.assertEqual(atom.resname, "LYS")
            self.assertEqual(atom.chain, "A")
            self.assertAlmostEqual(atom.coords[0], i + 1)
            self.assertAlmostEqual(atom.coords[1], 10 + i + 1)
            self.assertAlmostEqual(atom.coords[2], 20 + i + 1)

    def test_read_pdb_multiple_models(self):
        with mk_pdb_models(3) as tmp_pdb:
            models = pdb.read_pdb(tmp_pdb.name)
        self.assertIsInstance(models, list)
        self.assertEqual(len(models), 3)
        for atoms in models:
            self.assertIsInstance(atoms, AtomCollection)
            self.assertTrue(len(atoms), 10)

    def test_read_pdb_as_dict_single_model(self):
        with mk_pdb_models(1) as tmp_pdb:
            models = pdb.read_pdb(tmp_pdb.name, as_dict=True)
        self.assertIsInstance(models, dict)
        self.assertEqual(list(models.keys()), ["1"])
        for atoms in models.values():
            self.assertIsInstance(atoms, AtomCollection)
            self.assertTrue(len(atoms), 10)

    def test_read_pdb_as_dict_multiple_models(self):
        with mk_pdb_models(3) as tmp_pdb:
            models = pdb.read_pdb(tmp_pdb.name, as_dict=True)
        self.assertIsInstance(models, dict)
        self.assertEqual(list(models.keys()), ["1", "2", "3"])
        for atoms in models.values():
            self.assertIsInstance(atoms, AtomCollection)
            self.assertTrue(len(atoms), 10)

    def test_read_pdb_as_dict_no_model(self):
        with mk_pdb_no_model() as tmp_pdb:
            with self.assertRaisesRegex(
                pdb.InvalidPDBFormatError,
                r"can't initialize dictionary without model identifier \(no MODEL found\)",
            ):
                pdb.read_pdb(tmp_pdb.name, as_dict=True)

    def test_read_pdb_single_model(self):
        with mk_pdb_10_atoms() as tmp_pdb:
            atoms = pdb.read_pdb(tmp_pdb.name)
            self.assertIsInstance(atoms, AtomCollection)
            self.assertEqual(len(atoms), 10)

    def test_read_pdb_single_model_as_dict(self):
        with mk_pdb_10_atoms() as tmp_pdb:
            models = pdb.read_pdb(tmp_pdb.name, as_dict=True)
            self.assertIsInstance(models, dict)
            keys = list(models.keys())
            self.assertEqual(keys, ["1"])
            self.assertEqual(len(models["1"]), 10)
