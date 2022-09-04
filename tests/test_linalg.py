# Python core libraries.
import math

# Unit-test libraries.
from pytest import approx

from hypothesis import given
from hypothesis.strategies import composite, floats, integers
from hypothesis.extra.numpy import arrays

# Type-hinting special types.
from typing import Callable

# Scientific libraries.
import numpy as np

# PTools imports.
import ptools.linalg as linalg


def test_angle():
    u = (2, 2, 0)
    v = (0, 3, 0)
    angle = linalg.angle(u, v)
    assert math.degrees(angle) == approx(45)


@composite
def np_friendly_arrays(draw: Callable) -> np.array:
    shape = (
        draw(integers(min_value=2, max_value=8192)),
        draw(integers(min_value=2, max_value=8192)),
    )
    return draw(
        arrays(
            dtype=float, shape=shape, elements=floats(min_value=1e-10, max_value=1e10)
        )
    )


def test_centroid():
    x = np.random.rand(1000, 3)
    np.testing.assert_array_almost_equal(linalg.centroid(x), np.mean(x, axis=0))


def test_center_of_mass():
    def center_of_mass(x: np.ndarray, w: np.ndarray) -> np.ndarray:
        """Numpy-style COM calculation... way slower than PTools."""
        return np.average(x, axis=0, weights=w)

    x = np.random.rand(1000, 3)
    w = np.ones(x.shape[0])
    actual = linalg.center_of_mass(x, weights=w)
    expected = center_of_mass(x, w)
    np.testing.assert_array_almost_equal(actual, expected)




# def test_tensor_of_inertia_invalid_method():
#     err = (
#         r"parameter 'method' should be 'accurate' or 'fast' "
#         r"\(found method='foobar'\)"
#     )
#     with pytest.raises(ValueError, match=err):
#         linalg.tensor_of_inertia(np.zeros((3, 3)), method="foobar")


# def test_tensor_of_inertia_accurate_fails_without_weights(self):
#     array = np.array(((0, 0, 0), (1, 1, 1)), dtype=float)
#     err = "need weights to compute accurate tensor of inertia"
#     with self.assertRaisesRegex(ValueError, err):
#         spatial.tensor_of_inertia(array, method="accurate")

