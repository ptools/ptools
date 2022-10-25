"""Unit tests for ptools.coordinates."""

import random
import numpy as np

import pytest

from .testing.moreassert import assert_array_almost_equal

from ptools.coordinates import Coordinates3D, Invalid3DCoordinates


def generate_random_coordinates(shape: tuple[int, int] = (5, 3)) -> np.ndarray:
    """Generates a random np.ndarray."""
    return np.random.randn(*shape) * 100


def test_initialization_success():
    input_array = generate_random_coordinates()
    expected = input_array
    actual = Coordinates3D(input_array).array
    assert_array_almost_equal(expected, actual)



def test_initialization_from_vector_success():
    input_array = np.zeros(3)
    expected = input_array
    actual = Coordinates3D(input_array).array
    assert_array_almost_equal(expected, actual)



def test_initialization_failure():
    for ncol in (1, 2, 4, 5, 6):
        input_array = generate_random_coordinates(shape=(5, ncol))
        err = rf"cannot initialize 3D-coordinates from array with shape \(5, {ncol}\)"
        with pytest.raises(Invalid3DCoordinates, match=err):
            Coordinates3D(input_array)


def test_initialization_from_vector_failure():
    # 3 is the only allowed size for a vector.
    for size in (1, 2, 4, 5, 6):
        err = rf"cannot initialize 3D-coordinates from array with shape \({size},\)"
        with pytest.raises(Invalid3DCoordinates, match=err):
            Coordinates3D(generate_random_coordinates(shape=(size,)))


def test_assignment_success():
    expected = generate_random_coordinates()
    initial = generate_random_coordinates()
    coordinates = Coordinates3D(initial)
    coordinates.array = expected  # assigns coordinates array to new value
    actual = coordinates.array
    assert_array_almost_equal(expected, actual)


# def test_assignment_failure():
#     actual = Coordinates3D(generate_random_coordinates())
#     err = r"cannot initialize 3D-coordinates from array with shape \(5, 5\)"
#     with pytest.raises(Invalid3DCoordinates, match=err):
#         actual.array = generate_random_coordinates(shape=(5, 5))  # bad shape


def test_asarray():
    expected = generate_random_coordinates()
    coordinates = Coordinates3D(expected)
    actual = np.asarray(coordinates)
    assert isinstance(actual, np.ndarray)
    assert_array_almost_equal(expected, actual)


def test_can_compare_coordinates_with_numpy_arrays():
    expected = generate_random_coordinates()
    coordinates = Coordinates3D(expected)
    assert_array_almost_equal(expected, coordinates)


def test_shape():
    input_array = generate_random_coordinates()
    coordinates = Coordinates3D(input_array)
    actual = input_array.shape
    expected = coordinates.shape
    assert expected == actual


def test_getitem_int():
    input_array = generate_random_coordinates()
    coordinates = Coordinates3D(input_array)
    i = random.randint(0, input_array.shape[0] - 1)
    expected = input_array[i]
    actual = coordinates[i]
    assert isinstance(actual, Coordinates3D)
    assert_array_almost_equal(expected, actual)


def test_getitem_slice():
    input_array = generate_random_coordinates(shape=(5, 3))
    coordinates = Coordinates3D(input_array)
    expected = input_array[0:3]
    actual = coordinates[0:3]
    assert isinstance(actual, Coordinates3D)
    assert_array_almost_equal(expected, actual)


def set_item_int():
    input_array = generate_random_coordinates()
    coordinates = Coordinates3D(input_array)
    expected = [1, 2, 3]
    coordinates[0] = expected
    actual = coordinates[0]
    assert_array_almost_equal(expected, actual)


def set_item_int_from_Coordinates3D():
    input_array = generate_random_coordinates()
    coordinates = Coordinates3D(input_array)
    expected = Coordinates3D([1, 2, 3])
    coordinates[0] = expected
    actual = coordinates[0]
    assert_array_almost_equal(expected, actual)


def set_item_slice_slice():
    expected = generate_random_coordinates(shape=(2, 3))
    input_array = generate_random_coordinates(shape=(5, 3))
    coordinates = Coordinates3D(input_array)
    coordinates[0:2] = expected
    actual = coordinates[0:2]
    assert_array_almost_equal(expected, actual)


def set_item_slice_from_Coordinates3D():
    expected = Coordinates3D(generate_random_coordinates(shape=(2, 3)))
    input_array = generate_random_coordinates(shape=(5, 3))
    coordinates = Coordinates3D(input_array)
    coordinates[0:2] = expected
    actual = coordinates[0:2]
    assert_array_almost_equal(expected, actual)


def test_set_item_int_fails():
    input_array = generate_random_coordinates()
    coordinates = Coordinates3D(input_array)
    expected = [1, 2, 3, 4]  # bad shape
    err = r"could not broadcast input array from shape \(4,\) into shape \(3,\)"
    with pytest.raises(ValueError, match=err):
        coordinates[0] = expected
