"""test_io - Tests for `ptools.io` module."""

import os
import random
import string
import tempfile
import unittest

from ptools import io
from ptools.particlecollection import ParticleCollection
from ptools.io import pdb

from .testing.io import (
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
            self.assertEqual(atom.residue_index, 1)
            self.assertEqual(atom.residue_name, "LYS")
            self.assertEqual(atom.chain, "A")
            self.assertAlmostEqual(atom.coordinates[0], i + 1)
            self.assertAlmostEqual(atom.coordinates[1], 10 + i + 1)
            self.assertAlmostEqual(atom.coordinates[2], 20 + i + 1)

    def test_read_pdb_multiple_models(self):
        with mk_pdb_models(3) as tmp_pdb:
            models = pdb.read_pdb(tmp_pdb.name)
        self.assertIsInstance(models, list)
        self.assertEqual(len(models), 3)
        for atoms in models:
            self.assertIsInstance(atoms, ParticleCollection)
            self.assertTrue(len(atoms), 10)

    def test_read_pdb_as_dict_single_model(self):
        with mk_pdb_models(1) as tmp_pdb:
            models = pdb.read_pdb(tmp_pdb.name, as_dict=True)
        self.assertIsInstance(models, dict)
        self.assertEqual(list(models.keys()), ["1"])
        for atoms in models.values():
            self.assertIsInstance(atoms, ParticleCollection)
            self.assertTrue(len(atoms), 10)

    def test_read_pdb_as_dict_multiple_models(self):
        with mk_pdb_models(3) as tmp_pdb:
            models = pdb.read_pdb(tmp_pdb.name, as_dict=True)
        self.assertIsInstance(models, dict)
        self.assertEqual(list(models.keys()), ["1", "2", "3"])
        for atoms in models.values():
            self.assertIsInstance(atoms, ParticleCollection)
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
            self.assertIsInstance(atoms, ParticleCollection)
            self.assertEqual(len(atoms), 10)

    def test_read_pdb_single_model_as_dict(self):
        with mk_pdb_10_atoms() as tmp_pdb:
            models = pdb.read_pdb(tmp_pdb.name, as_dict=True)
            self.assertIsInstance(models, dict)
            keys = list(models.keys())
            self.assertEqual(keys, ["1"])
            self.assertEqual(len(models["1"]), 10)


class TestFileExists(unittest.TestCase):
    def test_assert_file_actually_exists(self):
        with tempfile.NamedTemporaryFile() as f:
            try:
                io.assert_file_exists(f.name)
            except (FileNotFoundError, IsADirectoryError) as exc:
                self.fail(
                    f"exception of type {type(exc).__name__} "
                    f"raised unexpectedly: '{exc}'"
                )

    def test_assert_file_does_not_exists(self):
        basename = "".join(random.choices(string.ascii_letters, k=8))
        path = "/foo/bar/baz/" + basename
        with self.assertRaisesRegex(FileNotFoundError, path):
            io.assert_file_exists(path)

    def test_assert_file_is_directory(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            with self.assertRaisesRegex(IsADirectoryError, tmpdirname):
                io.assert_file_exists(tmpdirname)

    def test_backup_if_exists(self):
        # Create empty source file.
        source = "/tmp/test_ptools_backup.tmp"
        create_empty_file(source)

        # Backup and check backup exists.
        io.backup_if_exists(source)
        self.assertTrue(os.path.exists(source + ".1"))

        # Remove backup file.
        os.remove(source + ".1")

    def test_backup_if_exists_two_times(self):
        # Create empty source file.
        source = "/tmp/test_ptools_backup.tmp"
        create_empty_file(source)

        # Backup and check backup exists.
        io.backup_if_exists(source)

        # Do it again (need to create source file again, since it's been
        # renamed at last step).
        create_empty_file(source)
        io.backup_if_exists(source)

        # Check two backup files have been created.
        self.assertTrue(os.path.exists(source + ".1"))
        self.assertTrue(os.path.exists(source + ".2"))

        # Tear down: remove backup files.
        os.remove(source + ".1")
        os.remove(source + ".2")


def create_empty_file(path):
    """Creates an empty file located at `path`."""
    with open(path, "wt", encoding="utf-8"):
        pass
