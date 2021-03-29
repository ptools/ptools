
import random
import sys
import unittest


from ptools.rigidbody import RigidBody
import ptools.superpose as superpose
from ptools.superpose import Screw
from ptools.spatial import coord3d

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
        m = np.zeros((4, 4))
        s = superpose.mat_trans_2_screw(m)
        self.assertTrue(isinstance(s, Screw))

        self.assertAlmostEqual(s.angle, 1.57079632679)
        self.assertAlmostEqual(s.normtranslation, 0.0)
        assert_array_almost_equal(s.unit, coord3d(1.0, 0.0, 0.0))
        assert_array_almost_equal(s.point, coord3d(0.0, 0.0, 0.0))


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

    def test_copy(self):
        target = self.s.copy()
        self.assertAlmostEqual(target.angle, self.s.angle)
        self.assertAlmostEqual(target.normtranslation, self.s.normtranslation)
        assert_array_almost_equal(target.unit, self.s.unit)
        assert_array_almost_equal(target.point, self.s.point)

        # Make sure changing one does not change the other
        target.angle += 12
        self.assertAlmostEqual(target.angle, self.s.angle + 12)

        target.unit += 3
        assert_array_almost_equal(target.unit, self.s.unit + 3)



def random_float():
    """Return a random floatting point number in the range
    [-max_float; +max_float]."""
    max_float = sys.float_info.max
    return random.randrange(-max_float, +max_float)

