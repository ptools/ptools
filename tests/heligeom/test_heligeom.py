import math
import os
import unittest

import numpy as np

from ptools import RigidBody
from ptools.heligeom import heli_analyze, heli_construct
from ptools.spatial import coord3d

from ..testing.moreassert import assert_array_equal, assert_array_almost_equal

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TEST_1A74_PROT_RED = os.path.join(TEST_DATA_DIR, "1A74_prot.red")
TEST_2GLSA = os.path.join(TEST_DATA_DIR, "2GLS_A.pdb")
TEST_2GLSB = os.path.join(TEST_DATA_DIR, "2GLS_B.pdb")
TEST_REF_2GLSAB_N6 = os.path.join(TEST_DATA_DIR, "ref_2GLSAB-N6.pdb")
TEST_REF_2GLSAB_N3_Z = os.path.join(TEST_DATA_DIR, "ref_2GLSAB-N3-Zalign.pdb")


def move_rigidbody(rb, x=0, y=0, z=0):
    out = rb.copy()
    out.moveby([x, y, z])
    return out


class TestHeligeomSimple(unittest.TestCase):
    def setUp(self):
        self.mono1 = RigidBody.from_pdb(TEST_1A74_PROT_RED)
        self.dx = 15
        self.mono2 = move_rigidbody(self.mono1, x=self.dx)

    def test_analyze_x_translate(self):
        hp = heli_analyze(self.mono1, self.mono2)
        self.assertAlmostEqual(hp.angle, 0.0)
        self.assertAlmostEqual(hp.normtranslation, self.dx)
        self.assertAlmostEqual(hp.unit[0], 1.0)
        self.assertAlmostEqual(hp.unit[1], 0.0)
        self.assertAlmostEqual(hp.unit[2], 0.0)

    def test_analyze_x_translate_rotate(self):
        point = coord3d(0, 0, 0)
        axis = coord3d(1, 0, 0)
        angle = math.pi / 4
        self.mono2.ab_rotate(point, axis, angle)

        hp = heli_analyze(self.mono1, self.mono2)
        self.assertAlmostEqual(hp.angle, angle)
        self.assertAlmostEqual(hp.normtranslation, self.dx)
        self.assertAlmostEqual(hp.unit[0], 1.0)
        self.assertAlmostEqual(hp.unit[1], 0.0)
        self.assertAlmostEqual(hp.unit[2], 0.0)

    def test_heli_construct(self):
        """Non-regression test for heligeom.heli_construct."""
        hp = heli_analyze(self.mono1, self.mono2)
        result = heli_construct(self.mono1, hp, N=15)
        reference_file = os.path.join(
            TEST_DATA_DIR, "test_heli_construct_simple_result.npy"
        )
        reference = np.load(reference_file)
        assert_array_equal(result.coords, reference)

    # Ignores W0703: Catching too general exception Exception
    # pylint: disable=W0703
    def test_Z_true_implemented(self):
        """Tests that using heligeom.extend with Z=true does not raises an error."""
        hp = heli_analyze(self.mono1, self.mono2)  # N is random
        try:
            heli_construct(self.mono1, hp, N=15, Z=True)
        except Exception as error:
            self.fail(
                f"heli_construct with Z=True unexpectedly raised an exception: '{error}'"
            )


class TestHeligeom(unittest.TestCase):
    def setUp(self):
        self.mono1 = RigidBody.from_pdb(TEST_2GLSA)
        self.mono2 = RigidBody.from_pdb(TEST_2GLSB)
        self.ref = RigidBody.from_pdb(TEST_REF_2GLSAB_N6)
        self.n_monomers = 6

    def test_hp_data(self):
        """Test heligeom.heli_analyze results"""
        hp = heli_analyze(self.mono1, self.mono2)

        self.assertAlmostEqual(hp.angle, 1.04719867)
        assert_array_almost_equal(hp.point, [0.000436, -0.000296, 0])
        assert_array_almost_equal(hp.unit, [8.47123119e-07, -2.80109302e-06, 1])

    def test_heli_construct(self):
        """Tests that heligeom.heli_construct"""
        hp = heli_analyze(self.mono1, self.mono2)
        result = heli_construct(self.mono1, hp, N=self.n_monomers)

        self.assertEqual(result.topdb(), self.ref.topdb())

    def test_heli_construct_Zalign(self):
        """Tests that heligeom.heli_construct"""

        ref = RigidBody.from_pdb(TEST_REF_2GLSAB_N3_Z)

        hp = heli_analyze(self.mono1, self.mono2)
        result = heli_construct(self.mono1, hp, N=3, Z=True)

        self.assertEqual(result.topdb(), ref.topdb())


if __name__ == "__main__":
    unittest.main()
