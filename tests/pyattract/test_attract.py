
"""test_attract - Tests for `pyptools.pyattract.attract`."""

import sys
import unittest

import pyptools.pyattract.attract as attract

from . import TEST_RECEPTOR_RED, TEST_LIGAND_RED


class CaptureStderrTest(unittest.TestCase):
    """Base class to run unit tests that capture standard error."""

    @property
    def error(self):
        sys.stderr.seek(0)
        return sys.stderr.read().decode()


class TestAttractCommandLine(CaptureStderrTest):

    def test_receptor_not_provided(self):
        with self.assertRaisesRegex(SystemExit, '2'):
            attract.parse_command_line('-l bar'.split())
        self.assertIn('the following arguments are required: -r/--receptor',
                      self.error)

    def test_ligand_not_provided(self):
        with self.assertRaisesRegex(SystemExit, '2'):
            attract.parse_command_line('-r foo'.split())
        self.assertIn('the following arguments are required: -l/--ligand',
                      self.error)

    def test_minimal_command_line(self):
        attract.parse_command_line('-r foo -l bar'.split())

    def test_unrecognized_option(self):
        with self.assertRaisesRegex(SystemExit, '2'):
            attract.parse_command_line('-r foo -l bar --foobar'.split())
        self.assertIn('unrecognized arguments: --foobar', self.error)

    def test_single_option(self):
        args = attract.parse_command_line('-r foo -l bar --single'.split())
        self.assertTrue(args.single)
        args = attract.parse_command_line('-r foo -l bar -s'.split())
        self.assertTrue(args.single)

    def test_ref_option(self):
        args = attract.parse_command_line('-r foo -l bar --ref foo.pdb'.split())
        self.assertEqual(args.ref, 'foo.pdb')


class TestAttract(CaptureStderrTest):

    def test_receptor_not_found(self):
        args = ['-r', 'foo', '-l', TEST_LIGAND_RED]
        with self.assertRaisesRegex(FileNotFoundError, "No such file: 'foo'"):
            attract.main(args)

    def test_ligand_not_found(self):
        args = ['-r', TEST_RECEPTOR_RED, '-l', 'bar']
        with self.assertRaisesRegex(FileNotFoundError, "No such file: 'bar'"):
            attract.main(args)

    def test_run(self):
        args = ['-r', TEST_RECEPTOR_RED, '-l', TEST_LIGAND_RED]
        attract.main(args)
        raise NotImplementedError('TODO')
