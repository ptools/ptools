"""Unit tests for ptools.coordinates."""

import numpy as np
import pytest

from ptools.array3d import array3d, Invalid3DArrayError

from .testing.moreassert import assert_array_almost_equal


def generate_random_coordinates(shape: tuple[int, int] = (5, 3)) -> np.ndarray:
    """Generates a random np.ndarray."""
    return np.random.randn(*shape) * 100


def test_initialization_success():
    input_array = generate_random_coordinates()
    expected = input_array
    actual = array3d(input_array)
    assert_array_almost_equal(expected, actual)


def test_initialization_from_vector_success():
    input_array = np.zeros(3)
    expected = input_array
    actual = array3d(input_array)
    assert_array_almost_equal(expected, actual)


def test_initialization_failure():
    for ncol in (1, 2, 4, 5, 6):
        input_array = generate_random_coordinates(shape=(5, ncol))
        err = rf"cannot initialize 3D-coordinates from array with shape \(5, {ncol}\)"
        with pytest.raises(Invalid3DArrayError, match=err):
            array3d(input_array)


def test_initialization_from_vector_failure():
    # 3 is the only allowed size for a vector.
    for size in (1, 2, 4, 5, 6):
        err = rf"cannot initialize 3D-coordinates from array with shape \({size},\)"
        with pytest.raises(Invalid3DArrayError, match=err):
            array3d(generate_random_coordinates(shape=(size,)))


def test_asarray():
    expected = generate_random_coordinates()
    coordinates = array3d(expected)
    actual = np.asarray(coordinates)
    assert isinstance(actual, np.ndarray)
    assert_array_almost_equal(expected, actual)


def test_can_compare_coordinates_with_numpy_arrays():
    expected = generate_random_coordinates()
    coordinates = array3d(expected)
    assert_array_almost_equal(expected, coordinates)


def test_set_item_int_fails():
    input_array = generate_random_coordinates()
    coordinates = array3d(input_array)
    expected = [1, 2, 3, 4]  # bad shape
    err = r"could not broadcast input array from shape \(4,\) into shape \(3,\)"
    with pytest.raises(ValueError, match=err):
        coordinates[0] = expected


def test_zeros():
    expected = (0, 0, 0)
    actual = array3d.zeros()
    assert_array_almost_equal(expected, actual)

    expected = ((0, 0, 0), (0, 0, 0))
    actual = array3d.zeros((2, 3))
    assert_array_almost_equal(expected, actual)


def test_zeros_fails():
    err = rf"cannot initialize 3D-coordinates from array with shape \(2,\)"
    with pytest.raises(Invalid3DArrayError, match=err):
        array3d.zeros(2)

    err = rf"cannot initialize 3D-coordinates from array with shape \(5, 2\)"
    with pytest.raises(Invalid3DArrayError, match=err):
        array3d.zeros((5, 2))
