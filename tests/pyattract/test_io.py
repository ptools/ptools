
"""test_io - tests for `pyattract.io`."""

import unittest

import ptools.pyattract.io as io

from . import (TEST_ATTRACT_PARAMS, TEST_ATTRACT_PARAMS_WITH_LIGAND,
               TEST_DUM_RED_CONTENT, TEST_DUM_PDB_CONTENT,
               TEST_AMINON_CONTENT)
from ..testing.io import random_filename, mk_tmp_file, mk_empty_file


with open(TEST_ATTRACT_PARAMS, 'rt') as f:
    TEST_ATTRACT_PARAMETER_CONTENT = f.read()


class TestAttractIO(unittest.TestCase):

    def test_read_aminon(self):
        tmpfile = mk_tmp_file(content=TEST_AMINON_CONTENT)
        params = io.read_aminon(tmpfile.name)
        tmpfile.close()
        self.assertEqual(len(params), 5)
        self.assertTrue(all([len(p) == 2 for p in params]))

        # Check radii read ok.
        self.assertEqual(params[0][0], 2.000)
        self.assertEqual(params[1][0], 1.900)
        self.assertEqual(params[2][0], 1.950)
        self.assertEqual(params[3][0], 1.900)
        self.assertEqual(params[4][0], 1.900)

        # Check amplitudes read ok.
        self.assertEqual(params[0][1], 1.000)
        self.assertEqual(params[1][1], 1.000)
        self.assertEqual(params[2][1], 2.000)
        self.assertEqual(params[3][1], 0.600)
        self.assertEqual(params[4][1], 0.600)

    def test_read_forcefield_from_reduced(self):
        tmpfile = mk_tmp_file(content=TEST_DUM_RED_CONTENT)
        ff = io.read_forcefield_from_reduced(tmpfile.name)
        tmpfile.close()
        self.assertEqual('attract1', ff)

    def test_read_forcefield_from_reduced_no_header(self):
        tmpfile = mk_tmp_file(content=TEST_DUM_PDB_CONTENT)
        err = 'reduced PDB file first line must be a HEADER line'
        with self.assertRaisesRegex(IOError, err):
            io.read_forcefield_from_reduced(tmpfile.name)
        tmpfile.close()

    def test_read_forcefield_from_reduced_misformatted_header(self):
        # Remove first line from proper RED file content and replace it
        # with a line that only contains HEADER with nothing following.
        lines = TEST_DUM_RED_CONTENT.splitlines()[1:]
        lines.insert(0, 'HEADER    ')

        # Make sure the proper exeception is raised.
        tmpfile = mk_tmp_file(content='\n'.join(lines))
        err = '.* cannot read force field name from first line .*'
        with self.assertRaisesRegex(IOError, err):
            io.read_forcefield_from_reduced(tmpfile.name)
        tmpfile.close()

    def test_read_forcefield_from_reduced_bad_ff(self):
        # Replace force field name with a bad one i.e. not defined
        # in `ptools.pyattract.PYATTRACT_FORCEFIELDS`.
        content = TEST_DUM_RED_CONTENT.replace('ATTRACT1', 'FOO')
        # Make sure the proper exeception is raised.
        tmpfile = mk_tmp_file(content=content)
        err = '.* invalid force field name .*'
        with self.assertRaisesRegex(ValueError, err):
            io.read_forcefield_from_reduced(tmpfile.name)
        tmpfile.close()

    def test_check_ff_version_match(self):
        tmpfile_receptor = mk_tmp_file(TEST_DUM_RED_CONTENT)
        tmpfile_ligand = mk_tmp_file(TEST_DUM_RED_CONTENT)
        ff = io.check_ff_version_match(tmpfile_receptor.name, tmpfile_ligand.name)
        self.assertEqual(ff, 'attract1')
        tmpfile_receptor.close()
        tmpfile_ligand.close()

    def test_check_ff_version_match_fail(self):
        content = TEST_DUM_RED_CONTENT.replace('ATTRACT1', 'ATTRACT2')
        tmpfile_receptor = mk_tmp_file(TEST_DUM_RED_CONTENT)
        tmpfile_ligand = mk_tmp_file(content)
        err = 'receptor and ligand force field names do not match'
        with self.assertRaisesRegex(ValueError, err):
            io.check_ff_version_match(tmpfile_receptor.name, tmpfile_ligand.name)
        tmpfile_receptor.close()
        tmpfile_ligand.close()


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
