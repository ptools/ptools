"""test_attract - Tests for `ptools.pyattract.attract`."""

import sys
import unittest

import pytest

from ptools.cli.ptools_cli import main as ptools_cli

from . import TEST_RECEPTOR_RED, TEST_LIGAND_RED


class CaptureStderrTest(unittest.TestCase):
    """Base class to run unit tests that capture standard error."""

    @property
    def error(self):
        sys.stderr.seek(0)
        return sys.stderr.read().decode()


class TestAttract(CaptureStderrTest):
    def test_receptor_not_found(self):
        args = ["attract", "-r", "foo", "-l", TEST_LIGAND_RED]
        with self.assertRaises(FileNotFoundError):
            ptools_cli(args)

    def test_ligand_not_found(self):
        args = ["attract", "-r", TEST_RECEPTOR_RED, "-l", "bar"]
        with self.assertRaises(FileNotFoundError):
            ptools_cli(args)
