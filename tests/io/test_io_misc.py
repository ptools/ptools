"""Tests for miscelleanous I/O functions."""

import os
import random
import string
import tempfile
import unittest

import ptools.io as io

from ..testing.io import create_empty_file


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
