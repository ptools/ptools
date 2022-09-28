# Python core libraries.
import math
import random

# Unit-test libraries.
import pytest
from pytest import approx

# Scientific libraries.
import numpy as np

# PTools imports.
from ptools import linalg

# More test-specific imports.
from .testing import assert_array_almost_equal


def generate_random_array(
    low: float = 0.0, high: float = 1.0, shape: int | tuple[int] = 1
) -> np.ndarray:
    return np.random.uniform(low, high, shape)


def generate_random_coordinates() -> np.ndarray:
    """Generate a coordinates-like array.

    The output array's size is random as well.
    """
    N = random.randint(2, 3000)
    return generate_random_array(-1000, 1000, (N, 3))


def generate_empty_coordinates() -> np.ndarray:
    """Returns an empty array."""
    return np.zeros((0, 3))


def test_angle():
    u = (2, 2, 0)
    v = (0, 3, 0)
    angle = linalg.angle(u, v)
    assert math.degrees(angle) == approx(45)


class TestTranslationMatrix:
    """Namespace that holds unit-tests for `ptools.linalg.translation_matrix`."""

    def test_vector(self):
        x = generate_random_array(shape=(3,))
        actual = linalg.matrix.translation_matrix(x)
        expected = np.array(
            [
                [1, 0, 0, x[0]],
                [0, 1, 0, x[1]],
                [0, 0, 1, x[2]],
                [0, 0, 0, 1],
            ]
        )
        assert_array_almost_equal(actual, expected)

    def test_scalar(self):
        x = random.random()
        actual = linalg.matrix.translation_matrix(x)
        expected = np.array(
            [
                [1, 0, 0, x],
                [0, 1, 0, x],
                [0, 0, 1, x],
                [0, 0, 0, 1],
            ]
        )
        assert_array_almost_equal(actual, expected)


class TestRotationMatrix:
    """Namespace that holds unit-tests for `ptools.linalg.rotation_matrix`."""

    def test_rotation_x_by_90(self):
        expected = [
            [1.00, 0.00, 0.00],
            [0.00, 0.00, -1.00],
            [0.00, 1.00, 0.00],
        ]
        actual = linalg.rotation_matrix([90, 0, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_x_by_180(self):
        expected = [
            [1.00, 0.00, 0.00],
            [0.00, -1.00, -0.00],
            [0.00, 0.00, -1.00],
        ]
        actual = linalg.rotation_matrix([180, 0, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_x_by_10(self):
        expected = [
            [1.0000000, 0.000000, 0.000000],
            [0.0000000, 0.984808, -0.173648],
            [0.0000000, 0.173648, 0.984808],
        ]
        actual = linalg.rotation_matrix([10, 0, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_y_by_90(self):
        expected = [
            [0.00, 0.00, 1.00],
            [0.00, 1.00, 0.00],
            [-1.00, 0.00, 0.00],
        ]
        actual = linalg.rotation_matrix([0, 90, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_y_by_180(self):
        expected = [
            [-1.00, 0.00, 0.00],
            [0.00, 1.00, 0.00],
            [-0.00, 0.00, -1.00],
        ]
        actual = linalg.rotation_matrix([0, 180, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_y_by_10(self):
        expected = [
            [0.9848078, 0.000000, 0.173648],
            [0.0000000, 1.000000, 0.000000],
            [-0.1736482, 0.000000, 0.984808],
        ]
        actual = linalg.rotation_matrix([0, 10, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_z_by_90(self):
        expected = [
            [0.00, -1.00, 0.00],
            [1.00, 0.00, 0.00],
            [0.00, 0.00, 1.00],
        ]
        actual = linalg.rotation_matrix([0, 0, 90])
        assert_array_almost_equal(actual, expected)

    def test_rotation_z_by_180(self):
        expected = [
            [-1.00, -0.00, 0.00],
            [0.00, -1.00, 0.00],
            [0.00, 0.00, 1.00],
        ]
        actual = linalg.rotation_matrix([0, 0, 180])
        assert_array_almost_equal(actual, expected)

    def test_rotation_z_by_10(self):
        expected = [
            [0.9848078, -0.173648, 0.000000],
            [0.1736482, 0.984808, 0.000000],
            [0.0000000, 0.000000, 1.000000],
        ]
        actual = linalg.rotation_matrix([0, 0, 10])
        assert_array_almost_equal(actual, expected)

    def test_rotation_xyz(self):
        """Tests rotation order."""
        angles = [10, 50, 80]
        by_x = linalg.rotation_matrix([angles[0], 0, 0])
        by_y = linalg.rotation_matrix([0, angles[1], 0])
        by_z = linalg.rotation_matrix([0, 0, angles[2]])

        actual = linalg.rotation_matrix(angles)
        expected = by_z.dot(by_y).dot(by_x)
        assert_array_almost_equal(actual, expected)

    def test_rotation_zyx(self):
        """Tests rotation order."""
        angles = [10, 10, 10]
        by_x = linalg.rotation_matrix([angles[0], 0, 0])
        by_y = linalg.rotation_matrix([0, angles[1], 0])
        by_z = linalg.rotation_matrix([0, 0, angles[2]])

        actual = linalg.rotation_matrix(angles, sequence="zyx")
        expected = by_x.dot(by_y).dot(by_z)
        assert_array_almost_equal(actual, expected)


# ======================================================================================
class TestCentroid:
    """Namespace that holds unit-tests for ptools.linalg.centroid."""

    def test_random_coordinates(self):
        x = generate_random_coordinates()
        assert_array_almost_equal(linalg.centroid(x), np.mean(x, axis=0))

    def test_empty_coordinates(self):
        x = generate_empty_coordinates()
        with pytest.warns(RuntimeWarning, match="Mean of empty slice"):
            assert np.isnan(linalg.centroid(x)).all()


# ======================================================================================
class TestCenterOfMass:
    """Namespace that holds unit-tests for ptools.linalg.center_of_mass."""

    def test_random_coordinates(self):
        x = generate_random_coordinates()
        w = generate_random_array(shape=x.shape[0])
        actual = linalg.center_of_mass(x, w)
        expected = np.average(x, axis=0, weights=w)
        assert_array_almost_equal(actual, expected)

    def test_empty_coordinates(self):
        x = generate_empty_coordinates()
        w = np.ones(x.shape[0])
        err = "cannot compute center of mass of empty array"
        with pytest.raises(ZeroDivisionError, match=err):
            linalg.center_of_mass(x, w)

    def test_empty_weights(self):
        x = generate_random_coordinates()
        w = np.ones(0)
        err = "input array and weights should be the same size"
        with pytest.raises(ValueError, match=err):
            linalg.center_of_mass(x, w)


# ======================================================================================

# pylint: disable=R0903
class TestInertiaTensor:
    """Namespace that holds unit-tests for ptools.linalg.inertia_tensor."""

    def test_fails_if_array_is_not_n_by_3(self):
        x = np.array(((0.0, 0.0), (1.0, 1.0)))
        w = np.ones((x.shape[0]))
        err = "inertia tensor can only be calculated on a N x 3 array"
        with pytest.raises(ValueError, match=err):
            linalg.inertia_tensor(x, w)
