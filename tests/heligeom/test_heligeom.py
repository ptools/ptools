
import math
import os
import unittest

import ptools
from ptools import RigidBody
from ptools.heligeom import extend, heli_analyze, heli_construct
from ptools.spatial import coord3d

from ..testing.moreassert import assert_array_equal


TEST_1A74_PROT_RED = os.path.join(os.path.dirname(__file__), 'data', '1A74_prot.red')


class TestHeligeom(unittest.TestCase):
    def test_analyze_x_translate(self):
        mono1 = RigidBody(TEST_1A74_PROT_RED)
        mono2 = RigidBody(TEST_1A74_PROT_RED)
        delta_x = 15
        tr = coord3d(delta_x, 0, 0)
        mono2.translate(tr)
        hp = heli_analyze(mono1, mono2)
        self.assertAlmostEqual(hp.angle, 0.0)
        self.assertAlmostEqual(hp.normtranslation, delta_x)
        self.assertAlmostEqual(hp.unit[0], 1.0)
        self.assertAlmostEqual(hp.unit[1], 0.0)
        self.assertAlmostEqual(hp.unit[2], 0.0)

    def test_analyze_x_translate_rotate(self):
        mono1 = RigidBody(TEST_1A74_PROT_RED)
        mono2 = RigidBody(TEST_1A74_PROT_RED)
        delta_x = 15
        tr = coord3d(delta_x, 0, 0)
        mono2.translate(tr)

        point = coord3d(0, 0, 0)
        axis = coord3d(1, 0, 0)
        angle = math.pi / 4
        mono2.ab_rotate(point, axis, angle)

        hp = heli_analyze(mono1, mono2)
        self.assertAlmostEqual(hp.angle, angle)
        self.assertAlmostEqual(hp.normtranslation, delta_x)
        self.assertAlmostEqual(hp.unit[0], 1.0)
        self.assertAlmostEqual(hp.unit[1], 0.0)
        self.assertAlmostEqual(hp.unit[2], 0.0)

    def test_heli_construct(self):
        """Tests that heligeom.heli_construct does the same thing as heligeom.extend."""
        mono1 = RigidBody(TEST_1A74_PROT_RED)
        mono2 = RigidBody(TEST_1A74_PROT_RED)
        delta_x = 15
        tr = coord3d(delta_x, 0, 0)
        mono2.translate(tr)

        hp = heli_analyze(mono1, mono2)

        target = extend(hp, mono1, N=15)
        result = heli_construct(mono1, hp, N=15)
        assert_array_equal(result.coords, target.coords)

    def test_Z_true_not_implemented(self):
        """Tests that using heligeom.extend with Z=true raises an error."""
        mono1 = RigidBody(TEST_1A74_PROT_RED)
        mono2 = RigidBody(TEST_1A74_PROT_RED)
        delta_x = 15
        tr = coord3d(delta_x, 0, 0)
        mono2.translate(tr)
        hp = heli_analyze(mono1, mono2)
        target = extend(hp, mono1, N=15)

        with self.assertRaises(NotImplementedError):
            result = heli_construct(mono1, hp, N=15, Z=True)





if __name__ == "__main__":
    unittest.main()
