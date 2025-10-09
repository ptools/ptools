# Scientific libraries.
import numpy as np

# PTools imports.
from ptools.linalg import matrix

# More test-specific imports.
from ..testing import assert_array_almost_equal


def test_translation_matrix():
    direction = [1.0, 2.0, 3.0]
    actual = matrix.translation_matrix(direction)
    expected = np.array(
        [
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 2.0],
            [0.0, 0.0, 1.0, 3.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    assert_array_almost_equal(actual, expected)


def test_rotation_matrix():
    angles = [90.0, 0.0, 0.0]
    actual = matrix.rotation_matrix(angles, degrees=True)
    expected = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 1.0, 0.0],
        ]
    )
    assert_array_almost_equal(actual, expected)


def test_transformation_matrix():
    translation = [1.0, 2.0, 3.0]
    rotation = [90.0, 0.0, 0.0]
    actual = matrix.transformation_matrix(translation, rotation)
    expected = np.array(
        [
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0, 2.0],
            [0.0, 1.0, 0.0, 3.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    assert_array_almost_equal(actual, expected)


def test_rotation_matrix_around_axis():
    axis = [1.0, 0.0, 0.0]
    amount = 90.0
    center = [0.0, 0.0, 0.0]
    actual = matrix.rotation_matrix_around_axis(axis, amount, center, degrees=True)
    expected = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    assert_array_almost_equal(actual, expected)


def test_ab_rotation_matrix():
    a = [0.0, 0.0, 0.0]
    b = [1.0, 0.0, 0.0]
    amount = 90.0
    actual = matrix.ab_rotation_matrix(a, b, amount, degrees=True)
    expected = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    assert_array_almost_equal(actual, expected)
