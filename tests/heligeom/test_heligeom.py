
import math
import os
import unittest

import pyptools
from pyptools import RigidBody
from pyptools.heligeom import heli_analyze
from pyptools.spatial import coord3d


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


if __name__ == "__main__":
    unittest.main()
