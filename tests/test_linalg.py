# Python core libraries.
from enum import Enum
import math
import random

# Unit-test libraries.
import pytest
from pytest import approx

# Scientific libraries.
import numpy as np

# PTools imports.
import ptools.linalg as linalg

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
    """Namespace that holds unit-tests for ptools.linalg.centroid."""

    def test_vector(self):
        x = generate_random_array(shape=(3, ))
        actual = linalg.matrix.translation_matrix(x)
        expected = np.array([
            [1, 0, 0, x[0]],
            [0, 1, 0, x[1]],
            [0, 0, 1, x[2]],
            [0, 0, 0, 1],
        ])
        assert_array_almost_equal(actual, expected)

    def test_scalar(self):
        x = random.random()
        actual = linalg.matrix.translation_matrix(x)
        expected = np.array([
            [1, 0, 0, x],
            [0, 1, 0, x],
            [0, 0, 1, x],
            [0, 0, 0, 1],
        ])
        assert_array_almost_equal(actual, expected)


class TestCentroid:
    """Namespace that holds unit-tests for ptools.linalg.centroid."""

    def test_random_coordinates(self):
        x = generate_random_coordinates()
        assert_array_almost_equal(linalg.centroid(x), np.mean(x, axis=0))

    def test_empty_coordinates(self):
        x = generate_empty_coordinates()
        with pytest.warns(RuntimeWarning, match="Mean of empty slice"):
            assert np.isnan(linalg.centroid(x)).all()


class TestCenterOfMass:
    """Namespace that holds unit-tests for ptools.linalg.center_of_mass."""

    @staticmethod
    def test_random_coordinates():
        x = generate_random_coordinates()
        w = generate_random_array(shape=x.shape[0])
        actual = linalg.center_of_mass(x, w)
        expected = np.average(x, axis=0, weights=w)
        assert_array_almost_equal(actual, expected)

    @staticmethod
    def test_empty_coordinates():
        x = generate_empty_coordinates()
        w = np.ones(x.shape[0])
        err = "cannot compute center of mass of empty array"
        with pytest.raises(ZeroDivisionError, match=err):
            linalg.center_of_mass(x, w)

    @staticmethod
    def test_empty_weights():
        x = generate_random_coordinates()
        w = np.ones(0)
        err = "input array and weights should be the same size"
        with pytest.raises(ValueError, match=err):
            linalg.center_of_mass(x, w)


class TestInertiaTensor:
    """Namespace that holds unit-tests for ptools.linalg.inertia_tensor."""

    # @staticmethod
    # def test_tensor_of_inertia():
    #     x = np.array(((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)))
    #     w = np.ones((x.shape[0]))
    #     expected = np.full((3, 3), 0.5)
    #     assert_array_almost_equal(linalg.tensor_of_inertia(x, w), expected)

    @staticmethod
    def test_fails_if_array_is_not_n_by_3():
        x = np.array(((0.0, 0.0), (1.0, 1.0)))
        w = np.ones((x.shape[0]))
        err = "inertia tensor can only be calculated on a N x 3 array"
        with pytest.raises(ValueError, match=err):
            linalg.inertia_tensor(x, w)
