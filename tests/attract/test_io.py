"""test_io - tests for `pyattract.io`."""

import unittest

import ptools.io.readers.attract as io

from ..generators import (
    generate_empty_file,
    generate_random_filename,
    generate_tmp_file,
)
from . import (
    TEST_AMINON_CONTENT,
    TEST_ATTRACT_PARAMS,
    TEST_ATTRACT_PARAMS_WITH_LIGAND,
    TEST_DUM_PDB_CONTENT,
    TEST_DUM_RED_CONTENT,
)

with open(TEST_ATTRACT_PARAMS, "rt", encoding="utf-8") as f:
    TEST_ATTRACT_PARAMETER_CONTENT = f.read()


class TestAttractIO(unittest.TestCase):
    def test_read_aminon(self):
        with generate_tmp_file(content=TEST_AMINON_CONTENT) as tmpfile:
            params = io.read_aminon(tmpfile.name)

        self.assertEqual(len(params), 5)
        self.assertTrue(all(len(p) == 2 for p in params))

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
        with generate_tmp_file(content=TEST_DUM_RED_CONTENT) as tmpfile:
            ff = io.read_forcefield_from_reduced(tmpfile.name)
            self.assertEqual("attract1", ff)

    def test_read_forcefield_from_reduced_no_header(self):
        with generate_tmp_file(content=TEST_DUM_PDB_CONTENT) as tmpfile:
            err = "reduced PDB file first line must be a HEADER line"
            with self.assertRaisesRegex(IOError, err):
                io.read_forcefield_from_reduced(tmpfile.name)

    def test_read_forcefield_from_reduced_misformatted_header(self):
        # Remove first line from proper RED file content and replace it
        # with a line that only contains HEADER with nothing following.
        lines = TEST_DUM_RED_CONTENT.splitlines()[1:]
        lines.insert(0, "HEADER    ")

        # Make sure the proper exeception is raised.
        with generate_tmp_file(content="\n".join(lines)) as tmpfile:
            err = ".* cannot read force field name from first line .*"
            with self.assertRaisesRegex(IOError, err):
                io.read_forcefield_from_reduced(tmpfile.name)

    def test_check_ff_version_match(self):
        with generate_tmp_file(TEST_DUM_RED_CONTENT) as tmpfile_receptor:
            with generate_tmp_file(TEST_DUM_RED_CONTENT) as tmpfile_ligand:
                ff = io.check_ff_version_match(
                    tmpfile_receptor.name, tmpfile_ligand.name
                )
                self.assertEqual(ff, "attract1")

    def test_check_ff_version_match_fail(self):
        content = TEST_DUM_RED_CONTENT.replace("ATTRACT1", "ATTRACT2")
        with generate_tmp_file(TEST_DUM_RED_CONTENT) as tmpfile_receptor:
            with generate_tmp_file(content) as tmpfile_ligand:
                err = "receptor and ligand force field names do not match"
                with self.assertRaisesRegex(ValueError, err):
                    io.check_ff_version_match(
                        tmpfile_receptor.name, tmpfile_ligand.name
                    )


class TestReadAttractParameters(unittest.TestCase):
    def test_file_does_not_exist(self):
        filename = generate_random_filename()
        err = f"No such file or directory: '{filename}'"
        with self.assertRaisesRegex(IOError, err):
            io.read_attract_parameter(filename)

    def test_cannot_read_number_of_minimizations(self):
        with generate_empty_file() as tmpfile:
            filename = tmpfile.name
            err = "Cannot read number of minimizations from attract parameter file"
            with self.assertRaisesRegex(ValueError, err):
                io.read_attract_parameter(filename)

    def test_file_contains_only_one_line(self):
        content = TEST_ATTRACT_PARAMETER_CONTENT.splitlines()[0]
        with generate_tmp_file(content=content) as tmpfile:
            filename = tmpfile.name
            err = "Unexpectedly reached end of attract parameter file"
            with self.assertRaisesRegex(ValueError, err):
                io.read_attract_parameter(filename)

    def test_cannot_read_rstk(self):
        lines = TEST_ATTRACT_PARAMETER_CONTENT.splitlines()[:2]
        lines[-1] = " ".join(lines[-1].split()[:-1])  # remove rstk from example
        content = "\n".join(lines)
        with generate_tmp_file(content=content) as tmpfile:
            filename = tmpfile.name
            err = "Cannot read rstk from attract parameter file"
            with self.assertRaisesRegex(ValueError, err):
                io.read_attract_parameter(filename)

    def test_cannot_read_too_many_minimizations(self):
        lines = TEST_ATTRACT_PARAMETER_CONTENT.splitlines()[:2]
        content = "\n".join(lines)
        with generate_tmp_file(content=content) as tmpfile:
            filename = tmpfile.name
            err = "Cannot read minimizations from attract parameter file: expected 4, found 0"
            with self.assertRaisesRegex(ValueError, err):
                io.read_attract_parameter(filename)

    def test_misformatted_minimization(self):
        lines = TEST_ATTRACT_PARAMETER_CONTENT.splitlines()[:3]
        lines[2] = " ".join(lines[-1].split()[:2])
        content = "\n".join(lines)
        with generate_tmp_file(content=content) as tmpfile:
            filename = tmpfile.name
            err = (
                "Cannot read minimization line from attract parameter file: "
                "expected at least 3 values, found 2"
            )
            with self.assertRaisesRegex(ValueError, err):
                io.read_attract_parameter(filename)

    def test_catches_ligand_lines(self):
        parameters = io.read_attract_parameter(TEST_ATTRACT_PARAMS_WITH_LIGAND)
        self.assertEqual(parameters.lignames[0], "foo")
        self.assertEqual(parameters.lignames[1], "bar")
        self.assertEqual(parameters.lignames[2], "baz")

    def test_read_parameters_ok(self):
        parameters = io.read_attract_parameter(TEST_ATTRACT_PARAMS)
        self.assertEqual(parameters.nbminim, 4)
        self.assertEqual(parameters.lignames, [])
        self.assertEqual(len(parameters.minimlist), 4)
        self.assertEqual(
            parameters.minimlist[0],
            {"maxiter": 30, "squarecutoff": 3000.00, "rstk": 0.00050},
        )
        self.assertEqual(
            parameters.minimlist[1],
            {"maxiter": 50, "squarecutoff": 500.00, "rstk": 0.0},
        )
        self.assertEqual(
            parameters.minimlist[2],
            {"maxiter": 100, "squarecutoff": 50.00, "rstk": 0.0},
        )
        self.assertEqual(
            parameters.minimlist[3],
            {"maxiter": 500, "squarecutoff": 1150.00, "rstk": 0.0},
        )
