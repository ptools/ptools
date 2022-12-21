import math
import os
import unittest

import numpy as np

from ptools import RigidBody
from ptools.heligeom import heli_analyze, heli_construct, dist_axis, contact
from ptools.io import to_pdb
from ptools import transform


from ..testing.moreassert import assert_array_equal, assert_array_almost_equal

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TEST_1A74_PROT_RED = os.path.join(TEST_DATA_DIR, "1A74_prot.red")
TEST_2GLSA = os.path.join(TEST_DATA_DIR, "2GLS_A.pdb")
TEST_2GLSB = os.path.join(TEST_DATA_DIR, "2GLS_B.pdb")
TEST_REF_2GLSAB_N6 = os.path.join(TEST_DATA_DIR, "ref_2GLSAB-N6.pdb")
TEST_REF_2GLSAB_N3_Z = os.path.join(TEST_DATA_DIR, "ref_2GLSAB-N3-Zalign.pdb")


def move_rigidbody(rb, x=0, y=0, z=0):
    out = rb.copy()
    transform.moveby(out, [x, y, z])
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
        point = np.array((0, 0, 0))
        axis = np.array((1, 0, 0))
        angle = math.pi / 4
        transform.ab_rotate(self.mono2, point, axis, angle, degrees=False)

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
        assert_array_almost_equal(result.coordinates, reference)

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
        self.hp = heli_analyze(self.mono1, self.mono2)

    def test_hp_data(self):
        """Test heligeom.heli_analyze results"""
        self.assertAlmostEqual(self.hp.angle, 1.04719867)
        assert_array_almost_equal(self.hp.point, [0.000436, -0.000296, 0])
        assert_array_almost_equal(self.hp.unit, [8.47123119e-07, -2.80109302e-06, 1])

    def test_heli_construct(self):
        """Tests that heligeom.heli_construct"""
        result = heli_construct(self.mono1, self.hp, N=self.n_monomers)
        self.assertEqual(to_pdb(result), to_pdb(self.ref))

    def test_heli_construct_Zalign(self):
        """Tests that heligeom.heli_construct"""
        ref = RigidBody.from_pdb(TEST_REF_2GLSAB_N3_Z)
        result = heli_construct(self.mono1, self.hp, N=3, Z=True)
        self.assertEqual(to_pdb(result), to_pdb(ref))

    def test_dist_axis(self):
        """Tests for heligeom.distAxis"""
        dmin, dmax = dist_axis(self.mono1, self.hp)
        self.assertAlmostEqual(dmin, 12.1986158)
        self.assertAlmostEqual(dmax, 73.5932897)

    def test_contact(self):
        """Test for heligeom.contact"""
        residues_in_contacts = contact(self.mono1, self.mono2)
        assert set(((57, 338), (32, 193), (37, 205))).issubset(residues_in_contacts)


if __name__ == "__main__":
    unittest.main()
