"""test_io - Tests for `ptools.io` module."""

import os
import string
import tempfile
import random
import unittest

import ptools.io as io
import ptools.io.pdb as pdb

from . import TEST_PDB, TEST_PDB_3MODELS, TEST_PDB_ATOM_NAMES
from .testing.io import random_filename


class TestPDBIO(unittest.TestCase):
    def test_read_pdb_does_not_exist(self):
        filename = random_filename()
        with self.assertRaises(FileNotFoundError):
            pdb.read_pdb(filename)

    def test_is_atom_line(self):
        self.assertTrue(pdb.is_atom_line("ATOM  <rest_of_the_line>"))
        self.assertTrue(pdb.is_atom_line("HETATM<rest_of_the_line>"))
        self.assertFalse(pdb.is_atom_line("ATOMI"))
        self.assertFalse(pdb.is_atom_line("ANISOU"))

    def test_read_pdb(self):
        atoms = pdb.read_pdb(TEST_PDB)
        self.assertEqual(len(atoms), 10)
        for i in range(10):
            atom = atoms[i]
            self.assertEqual(atom.index, i + 1)
            self.assertEqual(atom.name, TEST_PDB_ATOM_NAMES[i])
            self.assertEqual(atom.resid, 1)
            self.assertEqual(atom.resname, "LYS")
            self.assertEqual(atom.chain, "A")
            self.assertAlmostEqual(atom.coords[0], i + 1)
            self.assertAlmostEqual(atom.coords[1], 10 + i + 1)
            self.assertAlmostEqual(atom.coords[2], 20 + i + 1)

    def test_read_pdb_multiple_models(self):
        models = pdb.read_pdb(TEST_PDB_3MODELS)
        self.assertEqual(len(models), 3)
        for atoms in models:
            self.assertTrue(len(atoms), 10)


class TestFileExists(unittest.TestCase):
    def test_assert_file_actually_exists(self):
        f = tempfile.NamedTemporaryFile()
        try:
            io.assert_file_exists(f.name)
        except Exception as exc:
            self.fail(
                f"exception of type {type(exc).__name__} "
                f"raised unexpectedly: '{exc}'"
            )
        f.close()

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
    with open(path, "wt"):
        pass
