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


class TestCentroid:
    """Namespace that holds unit-tests for ptools.linalg.centroid."""

    def test_random_coordinates(self):
        x = generate_random_coordinates()
        np.testing.assert_array_almost_equal(linalg.centroid(x), np.mean(x, axis=0))

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
        np.testing.assert_array_almost_equal(actual, expected)

    @staticmethod
    def test_empty_coordinates():
        x = generate_empty_coordinates()
        w = np.ones(x.shape[0])
        with pytest.raises(
            ZeroDivisionError, match="cannot compute center of mass of empty array"
        ):
            linalg.center_of_mass(x, w)

    @staticmethod
    def test_empty_weights():
        x = generate_random_coordinates()
        w = np.ones(0)
        with pytest.raises(
            ValueError, match="input array and weights should be the same size"
        ):
            linalg.center_of_mass(x, w)


class TestTensorOfInertia:
    """Namespace that holds unit-tests for ptools.linalg.tensor_of_inertia."""

    @staticmethod
    def test_invalid_method():
        x = generate_random_array()
        with pytest.raises(ValueError):
            linalg.tensor_of_inertia(x, method="accurate")  # method invalid type: str

    @staticmethod
    def test_accurate_fails_without_weights():
        x = generate_random_array()
        err = "need weights to compute accurate tensor of inertia"
        with pytest.raises(ValueError, match=err):
            linalg.tensor_of_inertia(x, method=linalg.Method.ACCURATE)
