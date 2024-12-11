import math
from pathlib import Path
import unittest

import numpy as np
from pytest import approx

from ptools import RigidBody
from ptools.heligeom import heli_analyze, heli_construct
from ptools.measure import contacts_by_residue, minmax_distance_to_axis
from ptools.io import to_pdb
from ptools import transform


from ..testing.moreassert import assert_array_equal, assert_array_almost_equal

TEST_DATA_DIR = Path(__file__).parent / "data"
TEST_1A74_PROT_RED = TEST_DATA_DIR / "1A74_prot.red"
TEST_2GLSA = TEST_DATA_DIR / "2GLS_A.pdb"
TEST_2GLSB = TEST_DATA_DIR / "2GLS_B.pdb"
TEST_REF_2GLSAB_N6 = TEST_DATA_DIR / "ref_2GLSAB-N6.pdb"
TEST_REF_COORDS_2GLSAB_N6 = TEST_DATA_DIR / "ref_2GLSAB-N6.npy"
TEST_REF_2GLSAB_N3_Z = TEST_DATA_DIR / "ref_2GLSAB-N3-Zalign.pdb"
TEST_REF_COORDS_2GLSAB_N3_Z = TEST_DATA_DIR / "ref_2GLSAB-N3-Zalign.npy"


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
        assert hp.angle == approx(0.0)
        assert hp.normtranslation == approx(self.dx)
        assert hp.unit[0] == approx(1.0)
        assert hp.unit[1] == approx(0.0)
        assert hp.unit[2] == approx(0.0)

    def test_analyze_x_translate_rotate(self):
        point = np.array((0, 0, 0))
        axis = np.array((1, 0, 0))
        angle = math.pi / 4
        transform.ab_rotate(self.mono2, point, axis, angle, degrees=False)

        hp = heli_analyze(self.mono1, self.mono2)

        assert hp.angle == approx(angle)
        assert hp.normtranslation == approx(self.dx)
        assert hp.unit[0] == approx(1.0)
        assert hp.unit[1] == approx(0.0)
        assert hp.unit[2] == approx(0.0)

    def test_heli_construct(self):
        """Non-regression test for heligeom.heli_construct."""
        hp = heli_analyze(self.mono1, self.mono2)
        result = heli_construct(self.mono1, hp, N=3)
        reference_file = TEST_DATA_DIR / "test_heli_construct_simple_result.npy"
        reference = np.load(reference_file)
        assert result.coordinates == approx(reference)


class TestHeligeom(unittest.TestCase):
    def setUp(self):
        self.mono1 = RigidBody.from_pdb(TEST_2GLSA)
        self.mono2 = RigidBody.from_pdb(TEST_2GLSB)
        self.ref = RigidBody.from_pdb(TEST_REF_2GLSAB_N6)
        self.n_monomers = 6
        self.hp = heli_analyze(self.mono1, self.mono2)

    def test_hp_data(self):
        """Test heligeom.heli_analyze results"""
        assert self.hp.angle == approx(1.04719867)
        assert self.hp.point == approx((0.000436, -0.000296, 0), abs=1e-6)
        assert self.hp.unit == approx((8.47123119e-07, -2.80109302e-06, 1))

    def test_heli_construct(self):
        """Tests heligeom.heli_construct"""
        result = heli_construct(self.mono1, self.hp, N=self.n_monomers)
        ref_coords = np.load(TEST_REF_COORDS_2GLSAB_N6)

        assert to_pdb(result) == to_pdb(self.ref)
        assert result.coordinates == approx(ref_coords)

    def test_heli_construct_Zalign(self):
        """Tests heligeom.heli_construct z-alignment"""
        from ptools.io import write_pdb as write_pdb
        mono1 = self.mono1.copy()
        mono2 = self.mono2.copy()
        write_pdb(mono1, "mono1.pdb")
        write_pdb(mono2, "mono2.pdb")
        # Test structures are have helix axis oriented along z-axis,
        # so rotate them to a different orientation
        transform.rotate_by(mono1, (0,90,0))
        transform.rotate_by(mono2, (0,90,0))
        write_pdb(mono1, "mono1y.pdb")
        write_pdb(mono2, "mono2y.pdb")
        hp = heli_analyze(mono1, mono2)
        print("\n", hp)
        for k,v in hp.__dict__.items():
            print(k, v)
        result = heli_construct(mono1, hp, N=3, Z=True)
        write_pdb(result, "resultz.pdb")

        ref = RigidBody.from_pdb(TEST_REF_2GLSAB_N3_Z)
        ref_coords = np.load(TEST_REF_COORDS_2GLSAB_N3_Z)
        maxdiff = np.max(np.abs(result.coordinates - ref_coords))
        print("Max calculated - reference coord difference is ", maxdiff)

        #assert to_pdb(result) == to_pdb(ref)
        assert maxdiff < 5e-4 

    def test_dist_axis(self):
        """Tests for heligeom.distAxis"""
        dmin, dmax = minmax_distance_to_axis(self.mono1, self.hp.unit, center=self.hp.point)
        assert dmin == approx(12.1986158)
        assert dmax == approx(73.5932897)

    def test_contact(self):
        """Test for heligeom.contact"""
        residues_in_contacts = contacts_by_residue(self.mono1, self.mono2, 5)
        expected = set(((57, 338), (32, 193), (37, 205)))
        assert expected.issubset(residues_in_contacts)


if __name__ == "__main__":
    unittest.main()
