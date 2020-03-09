
import random
import sys
import unittest


from pyptools.rigidbody import RigidBody
import pyptools.superpose as superpose
from pyptools.superpose import Screw
from pyptools.spatial import coord3d

from scipy.spatial.transform import Rotation


from . import TEST_LIGAND
from .testing import assert_array_almost_equal


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


class TestScrew(unittest.TestCase):

    def setUp(self):
        self.s = Screw()

    def test_get_set_angle(self):
        self.assertEqual(self.s.angle, 0)
        value = random_float()
        self.s.angle = value
        self.assertAlmostEqual(self.s.angle, value)

    def test_get_set_normtranslation(self):
        self.assertEqual(self.s.normtranslation, 0)
        value = random_float()
        self.s.normtranslation = value
        self.assertAlmostEqual(self.s.normtranslation, value)

    def test_get_set_unit_vector(self):
        assert_array_almost_equal(self.s.unit, coord3d())
        u = coord3d((random_float(), random_float(), random_float()))
        self.s.unit = u
        assert_array_almost_equal(self.s.unit, u)

    def test_get_set_point(self):
        assert_array_almost_equal(self.s.point, coord3d())
        u = coord3d((random_float(), random_float(), random_float()))
        self.s.point = u
        assert_array_almost_equal(self.s.point, u)


def random_float():
    """Return a random floatting point number in the range
    [-max_float; +max_float]."""
    max_float = sys.float_info.max
    return random.randrange(-max_float, +max_float)

