# Unit-test libraries.
import unittest
from pytest import approx

from hypothesis import given
from hypothesis.strategies import composite, floats
from hypothesis.extra.numpy import arrays

# Type-hinting special types.
from typing import Callable
from hypothesis.strategies import SearchStrategy

# Scientific libraries.
import numpy as np
from scipy.spatial.transform import Rotation

# PTools imports.
from ptools import superpose
from ptools.rigidbody import RigidBody
from ptools.spatial import coord3d, transformation_matrix
from ptools.superpose import Screw

# Test-specific imports.
from . import TEST_LIGAND


def test_screw_initialization():
    screw = Screw()
    assert screw.unit == approx(np.zeros(3))
    assert screw.point == approx(np.zeros(3))
    assert screw.normtranslation == approx(0)
    assert screw.angle == approx(0)


def test_screw_copy():
    source = Screw()
    target = source.copy()
    assert target.unit == approx(source.unit)
    assert target.point == approx(source.point)
    assert target.normtranslation == approx(source.normtranslation)
    assert target.angle == approx(source.angle)

    # Makes sure modifying one does not modify the other.
    target.angle += 12
    target.unit += 3
    assert target.angle == approx(source.angle + 12)
    assert target.unit == approx(source.unit + 3)


@composite
def vectors3(
    draw: Callable[[SearchStrategy[float]], float],
    min_value: float = 1e-8,
    max_value: float = 1e8,
) -> np.ndarray:
    """Returns a custom vector of 3 floats with custom boundaries."""
    return draw(
        arrays(
            dtype=float,
            shape=(3,),
            elements=floats(min_value=min_value, max_value=max_value),
        )
    )


class TestSuperpose(unittest.TestCase):
    def setUp(self):
        self.target = RigidBody.from_pdb(TEST_LIGAND)

    @given(vectors3(), vectors3())
    def test_fit(self, translation_vector: np.ndarray, rotation_vector: np.ndarray):
        mobile = self.target.copy()

        # Random translation and rotation.
        rotation_matrix = Rotation.from_euler("xyz", rotation_vector).as_matrix()
        mobile.moveby(translation_vector)
        mobile.rotate(rotation_matrix)

        superpose.fit(mobile, self.target)
        assert mobile.coords == approx(self.target.coords, rel=1e-4)

    def test_rmsd(self):
        # RMSD with copy should be 0.0
        mobile = self.target.copy()
        assert superpose.rmsd(mobile, self.target) == approx(0.0)

        # RMSD after translation of 10 units should be 10.
        mobile.translate((10, 0, 0))
        assert superpose.rmsd(mobile, self.target) == approx(10.0)

        # RMSD after translation and fit should be 0.0.
        assert superpose.rmsd(mobile, self.target, do_fit=True) == approx(0.0)

    def test_mat_trans_to_screw(self):
        """Those are non-regression tests."""
        # abs(1 + a - b - c) > EPSILON
        m = np.zeros((4, 4))
        s = superpose.mat_trans_2_screw(m)
        assert isinstance(s, Screw)

        assert s.angle == approx(1.57079632679)
        assert s.normtranslation == approx(0.0)
        assert s.unit == approx(coord3d(1.0, 0.0, 0.0))
        assert s.point == approx(coord3d(0.0, 0.0, 0.0))

        # abs(1 - a + b - c) > EPSILON
        s = superpose.mat_trans_2_screw(
            transformation_matrix(rotation=np.array([0, -34, 0]))
        )
        assert s.angle == approx(-0.59341194)
        assert s.normtranslation == approx(0.0)
        assert s.unit == approx(coord3d(0.0, 1.0, 0.0))
        assert s.point == approx(coord3d(0.0, 0.0, 0.0))

        # abs(1 - a - b + c) > EPSILON
        s = superpose.mat_trans_2_screw(
            transformation_matrix(rotation=np.array([0, 0, 90]))
        )
        assert s.angle == approx(1.5707963)
        assert s.normtranslation == approx(0.0)
        assert s.unit == approx(coord3d(0.0, 0.0, 1.0))
        assert s.point == approx(coord3d(0.0, 0.0, 0.0))

        # angle = 0
        t = transformation_matrix(
            rotation=np.zeros(3), translation=np.array([-3, 2, 1])
        )
        s = superpose.mat_trans_2_screw(t)
        assert s.angle == approx(0.0)
        assert s.normtranslation == approx(3.74165738)
        assert s.unit == approx(coord3d(-0.80178373, 0.53452248, 0.26726124))
        assert s.point == approx(coord3d(0.0, 0.0, 0.0))

        t = transformation_matrix(
            rotation=np.array([11, 19, 87]), translation=np.array([-1, 8, 11])
        )
        s = superpose.mat_trans_2_screw(t)
        assert s.angle == approx(1.58731230)
        assert s.normtranslation == approx(10.95636347)
        assert s.unit == approx(coord3d(0.25480885, 0.07588361, 0.9640094))
        assert s.point == approx(coord3d(0.0, 3.30358643, 21.22781297))
