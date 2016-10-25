
"""test_rigidbody - Tests for `pyptools.rigidbody module."""

import unittest

from pyptools.rigidbody import RigidBody, AttractRigidBody

from . import TEST_PDB, TEST_RED
from .testing.io import mk_tmp_file
from .testing.moreassert import assert_array_equal

# The content of a reduced PDB (RED).
REDFILE = {
    'content': """\
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
    'categories': [1, 27, 28, 1, 23, 1, 23, 1, 27, 28],
    'charges': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
}

REDFILE_NO_CATEGORY = {
    'content': """\
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

REDFILE_NO_CHARGE = {
    'content': """\
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


class TestRigidBody(unittest.TestCase):

    def test_constructor(self):
        rb = RigidBody(TEST_PDB)
        self.assertEqual(len(rb), 10)

    def test_copy_constructor(self):
        rb1 = RigidBody(TEST_PDB)
        rb2 = RigidBody(rb1)
        self.assertEqual(len(rb1), len(rb2))
        ref_coords = rb1.coords.copy()
        rb1.coords.fill(0)
        assert_array_equal(rb2.coords, ref_coords)


class TestAttractRigidBody(unittest.TestCase):

    def test_constructor(self):
        tmpfile = mk_tmp_file(content=REDFILE['content'])
        arb = AttractRigidBody(tmpfile.name)
        tmpfile.close()
        self.assertEqual(len(arb), 10)
        # !! atom categories store in AttractRigidBody are array indices,
        # !! therefore minored by one compared to what is in the RED file.
        assert_array_equal(arb.atom_categories + 1, REDFILE['categories'])
        assert_array_equal(arb.atom_charges, REDFILE['charges'])

    def test_constructor_fails_no_category(self):
        tmpfile = mk_tmp_file(content=REDFILE_NO_CATEGORY['content'])
        err = 'cannot initialize atom category'
        with self.assertRaisesRegex(IOError, err):
            AttractRigidBody(tmpfile.name)
        tmpfile.close()

    def test_constructor_fails_no_charge(self):
        tmpfile = mk_tmp_file(content=REDFILE_NO_CHARGE['content'])
        err = 'cannot initialize atom charge'
        with self.assertRaisesRegex(IOError, err):
            AttractRigidBody(tmpfile.name)
        tmpfile.close()

    def test_reset_forces(self):
        arb = AttractRigidBody(TEST_RED)
        assert_array_equal(arb.atom_forces, 0.0)
        arb.atom_forces[0] = [1, 2, 3]
        assert_array_equal(arb.atom_forces[0], [1, 2, 3])
        arb.reset_forces()
        assert_array_equal(arb.atom_forces, 0.0)

