
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
    def setUp(self):
        self.mono1 = RigidBody(TEST_1A74_PROT_RED)
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
        """Tests that heligeom.heli_construct does the same thing as heligeom.extend."""
        hp = heli_analyze(self.mono1, self.mono2)

        target = extend(hp, self.mono1, N=15)  # N is random
        result = heli_construct(self.mono1, hp, N=15)
        assert_array_equal(result.coords, target.coords)

    def test_Z_true_implemented(self):
        """Tests that using heligeom.extend with Z=true does not raises an error."""
        hp = heli_analyze(self.mono1, self.mono2)  # N is random
        target = extend(hp, self.mono1, N=15)
        try:
            result = heli_construct(self.mono1, hp, N=15, Z=True)
        except:
            self.fail("heli_construct with Z=True unexpectedly raised an exception")


def move_rigidbody(rb, x=0, y=0, z=0):
    out = rb.copy()
    out.moveby([x, y, z])
    return out


if __name__ == "__main__":
    unittest.main()
