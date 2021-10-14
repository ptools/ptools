"""test_rigidbody - Tests for `ptools.rigidbody module."""

import unittest

from ptools.rigidbody import RigidBody, AttractRigidBody
from ptools.io import InvalidPDBFormatError

from . import TEST_PDB, TEST_RED
from .testing.io import mk_tmp_file
from .testing.moreassert import assert_array_equal


# The content of a reduced PDB (RED).
REDFILE = {
    "content": """\
ATOM      1  CA  TYR     1       0.880   6.552  -1.114    1   0.000 0 0
ATOM      2  CSE TYR     1      -0.313   5.011  -0.781   27   0.000 0 0
ATOM      3  CSE TYR     1      -0.074   2.173  -0.857   28   0.000 0 0
ATOM      4  CA  SER     2       3.501   7.012   1.646    1   0.000 0 0
ATOM      5  CSE SER     2       3.640   8.268   1.516   23   0.000 0 0
ATOM      6  CA  SER     3       1.143   5.686   4.387    1   0.000 0 0
ATOM      7  CSE SER     3       0.555   6.008   5.467   23   0.000 0 0
ATOM      8  CA  TYR     4      -1.411   2.949   3.745    1   0.000 0 0
ATOM      9  CSE TYR     4      -1.675   0.958   3.911   27   0.000 0 0
ATOM     10  CSE TYR     4      -3.525  -1.079   3.005   28   0.000 0 0
""",
    "categories": [1, 27, 28, 1, 23, 1, 23, 1, 27, 28],
    "charges": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
}

REDFILE_NO_CATEGORY = {
    "content": """\
ATOM      1  CA  TYR     1       0.880   6.552  -1.114
ATOM      2  CSE TYR     1      -0.313   5.011  -0.781
ATOM      3  CSE TYR     1      -0.074   2.173  -0.857
ATOM      4  CA  SER     2       3.501   7.012   1.646
ATOM      5  CSE SER     2       3.640   8.268   1.516
ATOM      6  CA  SER     3       1.143   5.686   4.387
ATOM      7  CSE SER     3       0.555   6.008   5.467
ATOM      8  CA  TYR     4      -1.411   2.949   3.745
ATOM      9  CSE TYR     4      -1.675   0.958   3.911
ATOM     10  CSE TYR     4      -3.525  -1.079   3.005
""",
}


REDFILE_NO_CHARGE = {
    "content": """\
ATOM      1  CA  TYR     1       0.880   6.552  -1.114    1
ATOM      2  CSE TYR     1      -0.313   5.011  -0.781   27
ATOM      3  CSE TYR     1      -0.074   2.173  -0.857   28
ATOM      4  CA  SER     2       3.501   7.012   1.646    1
ATOM      5  CSE SER     2       3.640   8.268   1.516   23
ATOM      6  CA  SER     3       1.143   5.686   4.387    1
ATOM      7  CSE SER     3       0.555   6.008   5.467   23
ATOM      8  CA  TYR     4      -1.411   2.949   3.745    1
ATOM      9  CSE TYR     4      -1.675   0.958   3.911   27
ATOM     10  CSE TYR     4      -3.525  -1.079   3.005   28
""",
}


REDFILE_INVALID_CATEGORY = {
    "content": """\
ATOM      1  CA  TYR     1       0.880   6.552  -1.114  1.2   0.000 0 0
ATOM      2  CSE TYR     1      -0.313   5.011  -0.781 27.2   0.000 0 0
ATOM      3  CSE TYR     1      -0.074   2.173  -0.857 28.2   0.000 0 0
ATOM      4  CA  SER     2       3.501   7.012   1.646  1.2   0.000 0 0
ATOM      5  CSE SER     2       3.640   8.268   1.516 23.2   0.000 0 0
ATOM      6  CA  SER     3       1.143   5.686   4.387  1.2   0.000 0 0
ATOM      7  CSE SER     3       0.555   6.008   5.467 23.2   0.000 0 0
ATOM      8  CA  TYR     4      -1.411   2.949   3.745  1.2   0.000 0 0
ATOM      9  CSE TYR     4      -1.675   0.958   3.911 27.2   0.000 0 0
ATOM     10  CSE TYR     4      -3.525  -1.079   3.005 28.2   0.000 0 0
""",
}


REDFILE_INVALID_CHARGES = {
    "content": """\
ATOM      1  CA  TYR     1       0.880   6.552  -1.114    1   TYR   0 0
ATOM      2  CSE TYR     1      -0.313   5.011  -0.781   27   TYR   0 0
ATOM      3  CSE TYR     1      -0.074   2.173  -0.857   28   TYR   0 0
ATOM      4  CA  SER     2       3.501   7.012   1.646    1   SER   0 0
ATOM      5  CSE SER     2       3.640   8.268   1.516   23   SER   0 0
ATOM      6  CA  SER     3       1.143   5.686   4.387    1   SER   0 0
ATOM      7  CSE SER     3       0.555   6.008   5.467   23   SER   0 0
ATOM      8  CA  TYR     4      -1.411   2.949   3.745    1   TYR   0 0
ATOM      9  CSE TYR     4      -1.675   0.958   3.911   27   TYR   0 0
ATOM     10  CSE TYR     4      -3.525  -1.079   3.005   28   TYR   0 0
""",
}




class TestRigidBody(unittest.TestCase):
    def setUp(self):
        self.rb = RigidBody(TEST_PDB)

    def test_constructor(self):
        self.assertEqual(len(self.rb), 10)

    def _assert_copy_successful(self, thecopy):
        # Check that both RigidBody instances have the same number of atoms.
        self.assertEqual(len(thecopy), len(self.rb))
        self.assertEqual(thecopy.coords.shape, (10, 3))

        # Change parent AtomCollection coordinates and make sure it does not
        # affect the copy.
        ref_coords = self.rb.coords.copy()
        self.rb.coords.fill(0)
        assert_array_equal(thecopy.coords, ref_coords)

    def test_copy_constructor1(self):
        thecopy = self.rb.copy()
        self._assert_copy_successful(thecopy)

    def test_copy_constructor2(self):
        thecopy = RigidBody.copy(self.rb)
        self._assert_copy_successful(thecopy)

    def test_getitem_returns_rigidbody(self):
        atoms = self.rb[:5]
        self.assertIsInstance(atoms, RigidBody)


class TestAttractRigidBody(unittest.TestCase):
    def test_constructor(self):
        with mk_tmp_file(content=REDFILE["content"]) as tmpfile:
            arb = AttractRigidBody(tmpfile.name)

        self.assertEqual(len(arb), 10)
        # !! atom categories store in AttractRigidBody are array indices,
        # !! therefore minored by one compared to what is in the RED file.
        assert_array_equal(arb.atom_categories + 1, REDFILE["categories"])
        assert_array_equal(arb.atom_charges, REDFILE["charges"])

    def test_constructor_fails_no_categories(self):
        with mk_tmp_file(content=REDFILE_NO_CATEGORY["content"]) as tmpfile:
            err = "Expected atom categories and charges, found"
            with self.assertRaisesRegex(InvalidPDBFormatError, err):
                AttractRigidBody(tmpfile.name)

    def test_constructor_fails_no_charges(self):
        with mk_tmp_file(content=REDFILE_NO_CATEGORY["content"]) as tmpfile:
            err = "Expected atom categories and charges, found"
            with self.assertRaisesRegex(InvalidPDBFormatError, err):
                AttractRigidBody(tmpfile.name)

    def test_constructor_fails_invalid_categories(self):
        with mk_tmp_file(content=REDFILE_INVALID_CATEGORY["content"]) as tmpfile:
            err = "Atom category expects an int"
            with self.assertRaisesRegex(InvalidPDBFormatError, err):
                AttractRigidBody(tmpfile.name)

    def test_constructor_fails_invalid_charges(self):
        with mk_tmp_file(content=REDFILE_INVALID_CHARGES["content"]) as tmpfile:
            err = "Atom charge expects a float"
            with self.assertRaisesRegex(InvalidPDBFormatError, err):
                AttractRigidBody(tmpfile.name)


    # Ignores R0201: Method could be a function (no-self-use)
    # pylint: disable=R0201
    def test_reset_forces(self):
        arb = AttractRigidBody(TEST_RED)
        assert_array_equal(arb.atom_forces, 0.0)
        arb.atom_forces[0] = [1, 2, 3]
        assert_array_equal(arb.atom_forces[0], [1, 2, 3])
        arb.reset_forces()
        assert_array_equal(arb.atom_forces, 0.0)
