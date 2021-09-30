import random
import sys
import unittest

import numpy as np
from scipy.spatial.transform import Rotation

from ptools import superpose
from ptools.rigidbody import RigidBody
from ptools.spatial import coord3d, transformation_matrix
from ptools.superpose import Screw

from . import TEST_LIGAND
from .testing import assert_array_almost_equal


class TestScrew(unittest.TestCase):
    def setUp(self):
        self.screw = Screw()

    def test_get_set_angle(self):
        self.assertEqual(self.screw.angle, 0)
        value = random_float()
        self.screw.angle = value
        self.assertAlmostEqual(self.screw.angle, value)

    def test_get_set_normtranslation(self):
        self.assertEqual(self.screw.normtranslation, 0)
        value = random_float()
        self.screw.normtranslation = value
        self.assertAlmostEqual(self.screw.normtranslation, value)

    def test_get_set_unit_vector(self):
        assert_array_almost_equal(self.screw.unit, coord3d())
        u = coord3d((random_float(), random_float(), random_float()))
        self.screw.unit = u
        assert_array_almost_equal(self.screw.unit, u)

    def test_get_set_point(self):
        assert_array_almost_equal(self.screw.point, coord3d())
        u = coord3d((random_float(), random_float(), random_float()))
        self.screw.point = u
        assert_array_almost_equal(self.screw.point, u)

    def test_copy(self):
        target = self.screw.copy()
        self.assertAlmostEqual(target.angle, self.screw.angle)
        self.assertAlmostEqual(target.normtranslation, self.screw.normtranslation)
        assert_array_almost_equal(target.unit, self.screw.unit)
        assert_array_almost_equal(target.point, self.screw.point)

        # Make sure changing one does not change the other
        target.angle += 12
        self.assertAlmostEqual(target.angle, self.screw.angle + 12)

        target.unit += 3
        assert_array_almost_equal(target.unit, self.screw.unit + 3)


class TestSuperpose(unittest.TestCase):
    def setUp(self):
        self.target = RigidBody(TEST_LIGAND)

    def test_fit(self):
        mobile = self.target.copy()

        # Random translation and rotation.
        t = [-4.2, 5.1, 20.2]
        r = Rotation.from_euler("xyz", [-10, 30, 90]).as_matrix()

        mobile.moveby(t)
        mobile.rotate(r)

        superpose.fit(mobile, self.target)

        assert_array_almost_equal(mobile.coords, self.target.coords)

    def test_rmsd(self):
        # RMSD with copy should be 0.0
        mobile = self.target.copy()
        self.assertAlmostEqual(superpose.rmsd(mobile, self.target), 0.0)

        # RMSD after translation of 10 units should be 10.
        mobile.translate((10, 0, 0))
        self.assertAlmostEqual(superpose.rmsd(mobile, self.target), 10.0)

        # RMSD after translation and fit should be 0.0.
        self.assertAlmostEqual(superpose.rmsd(mobile, self.target, do_fit=True), 0)

    def test_mat_trans_to_screw(self):

        # abs(1 + a - b - c) > EPSILON
        m = np.zeros((4, 4))
        s = superpose.mat_trans_2_screw(m)
        self.assertTrue(isinstance(s, Screw))

        self.assertAlmostEqual(s.angle, 1.57079632679)
        self.assertAlmostEqual(s.normtranslation, 0.0)
        assert_array_almost_equal(s.unit, coord3d(1.0, 0.0, 0.0))
        assert_array_almost_equal(s.point, coord3d(0.0, 0.0, 0.0))

        # abs(1 - a + b - c) > EPSILON
        s = superpose.mat_trans_2_screw(
            transformation_matrix(rotation=np.array([0, -34, 0]))
        )
        self.assertAlmostEqual(s.angle, -0.59341194)
        self.assertAlmostEqual(s.normtranslation, 0.0)
        assert_array_almost_equal(s.unit, coord3d(0.0, 1.0, 0.0))
        assert_array_almost_equal(s.point, coord3d(0.0, 0.0, 0.0))

        # abs(1 - a - b + c) > EPSILON
        s = superpose.mat_trans_2_screw(
            transformation_matrix(rotation=np.array([0, 0, 90]))
        )
        self.assertAlmostEqual(s.angle, 1.5707963)
        self.assertAlmostEqual(s.normtranslation, 0.0)
        assert_array_almost_equal(s.unit, coord3d(0.0, 0.0, 1.0))
        assert_array_almost_equal(s.point, coord3d(0.0, 0.0, 0.0))

        # angle = 0
        t = transformation_matrix(
            rotation=np.zeros(3), translation=np.array([-3, 2, 1])
        )
        s = superpose.mat_trans_2_screw(t)
        self.assertAlmostEqual(s.angle, 0.0)
        self.assertAlmostEqual(s.normtranslation, 3.74165738)
        assert_array_almost_equal(s.unit, coord3d(-0.80178373, 0.53452248, 0.26726124))
        assert_array_almost_equal(s.point, coord3d(0.0, 0.0, 0.0))

        t = transformation_matrix(
            rotation=np.array([11, 19, 87]), translation=np.array([-1, 8, 11])
        )
        s = superpose.mat_trans_2_screw(t)
        self.assertAlmostEqual(s.angle, 1.58731230)
        self.assertAlmostEqual(s.normtranslation, 10.95636347)
        assert_array_almost_equal(s.unit, coord3d(0.25480885, 0.07588361, 0.9640094))
        assert_array_almost_equal(s.point, coord3d(0.0, 3.30358643, 21.22781297))


def random_float():
    """Return a random floatting point number in the range
    [-max_float; +max_float]."""
    max_float = sys.float_info.max
    return random.randrange(-max_float, +max_float)
