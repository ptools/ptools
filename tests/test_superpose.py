
import random
import unittest

from pyptools.rigidbody import RigidBody
import pyptools.superpose as superpose

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



