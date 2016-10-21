
"""test_io - tests for `pyattract.io`."""

import unittest

import pyptools.pyattract.io as io

from . import TEST_ATTRACT_PARAMS, TEST_ATTRACT_PARAMS_WITH_LIGAND
from ..testing.io import random_filename, mk_tmp_file, mk_empty_file


with open(TEST_ATTRACT_PARAMS, 'rt') as f:
    TEST_ATTRACT_PARAMETER_CONTENT = f.read()


class TestReadAttractParameters(unittest.TestCase):

    def test_file_does_not_exist(self):
        filename = random_filename()
        err = "No such file or directory: '{}'".format(filename)
        with self.assertRaisesRegex(IOError, err):
            io.read_attract_parameter(filename)

    def test_cannot_read_number_of_minimizations(self):
        tmpfile = mk_empty_file()
        filename = tmpfile.name
        err = 'Cannot read number of minimizations from attract parameter file'
        with self.assertRaisesRegex(ValueError, err):
            io.read_attract_parameter(filename)
        tmpfile.close()

    def test_file_contains_only_one_line(self):
        content = TEST_ATTRACT_PARAMETER_CONTENT.splitlines()[0]
        tmpfile = mk_tmp_file(content=content)
        filename = tmpfile.name
        err = 'Unexpectedly reached end of attract parameter file'
        with self.assertRaisesRegex(ValueError, err):
            io.read_attract_parameter(filename)
        tmpfile.close()

    def test_cannot_read_rstk(self):
        lines = TEST_ATTRACT_PARAMETER_CONTENT.splitlines()[:2]
        lines[-1] = ' '.join(lines[-1].split()[:-1])  # remove rstk from example
        content = '\n'.join(lines)
        tmpfile = mk_tmp_file(content=content)
        filename = tmpfile.name
        err = 'Cannot read rstk from attract parameter file'
        with self.assertRaisesRegex(ValueError, err):
            io.read_attract_parameter(filename)
        tmpfile.close()

    def test_cannot_read_too_many_minimizations(self):
        lines = TEST_ATTRACT_PARAMETER_CONTENT.splitlines()[:2]
        content = '\n'.join(lines)
        tmpfile = mk_tmp_file(content=content)
        filename = tmpfile.name
        err = 'Cannot read minimizations from attract parameter file: expected 4, found 0'
        with self.assertRaisesRegex(ValueError, err):
            io.read_attract_parameter(filename)
        tmpfile.close()

    def test_misformatted_minimization(self):
        lines = TEST_ATTRACT_PARAMETER_CONTENT.splitlines()[:3]
        lines[2] = ' '.join(lines[-1].split()[:2])
        content = '\n'.join(lines)
        tmpfile = mk_tmp_file(content=content)
        filename = tmpfile.name
        err = 'Cannot read minimization line from attract parameter file: '\
              'expected at least 3 values, found 2'
        with self.assertRaisesRegex(ValueError, err):
            io.read_attract_parameter(filename)
        tmpfile.close()

    def test_catches_ligand_lines(self):
        parameters = io.read_attract_parameter(TEST_ATTRACT_PARAMS_WITH_LIGAND)
        lignames = parameters[1]
        self.assertEqual(lignames[0], 'foo')
        self.assertEqual(lignames[1], 'bar')
        self.assertEqual(lignames[2], 'baz')

    def test_read_parameters_ok(self):
        nbminim, lignames, minimlist, rstk = io.read_attract_parameter(TEST_ATTRACT_PARAMS)
        self.assertEqual(nbminim, 4)
        self.assertEqual(lignames, [])
        self.assertEqual(len(minimlist), 4)
        self.assertEqual(minimlist[0], {'maxiter': 30, 'squarecutoff': 3000.00, 'rstk': 0.00050})
        self.assertEqual(minimlist[1], {'maxiter': 50, 'squarecutoff': 500.00, 'rstk': 0.0})
        self.assertEqual(minimlist[2], {'maxiter': 100, 'squarecutoff': 50.00, 'rstk': 0.0})
        self.assertEqual(minimlist[3], {'maxiter': 500, 'squarecutoff': 1150.00, 'rstk': 0.0})
