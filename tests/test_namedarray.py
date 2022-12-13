""""test_atomproperty - Tests for ``NamedArray``."""

from ptools.namedarray import NamedArray

import numpy as np
import pytest

from .testing import assert_array_almost_equal, assert_array_not_almost_equal


def test_initialization_converts_to_numpy():
    prop = NamedArray("foo", "bar", (1, 2, 3))
    assert isinstance(prop.values, np.ndarray)


def test_assignment_converts_to_numpy():
    prop = NamedArray("foo", "bar", np.ones(5))
    prop.values = (1, 2, 3)
    assert isinstance(prop.values, np.ndarray)

    prop.values = 1
    assert isinstance(prop.values, np.ndarray)


def test_assignment_does_actual_copies():
    prop = NamedArray("foo", "bar", np.ones(5))
    expected = np.array((1, 2, 3))

    prop.values = expected
    assert_array_almost_equal(prop.values, expected)

    expected[0] = 42
    expected[1] = 17
    expected[2] = 256

    assert_array_not_almost_equal(prop.values, expected)


def test_equality():
    left = NamedArray("foo", "bar", np.ones(5))
    right = NamedArray("foo", "bar", np.ones(5))
    assert left == right

    right = NamedArray("foo", "bar", np.zeros(5))
    assert left != right


def test_equality_float():
    left = NamedArray("foo", "bar", np.ones(5))
    right = NamedArray("foo", "bar", np.ones(5))
    assert left == right

    right = NamedArray("foo", "bar", np.ones(5) + 1e-6)
    assert left == right

    right = NamedArray("foo", "bar", np.ones(5) + 1e-4)
    assert left != right

    # Check array comparison function works as expected
    left._array_comparison_func = lambda x, y: True
    assert left == right


def test_equality_type_error():
    left = NamedArray("foo", "bar", np.ones(5))
    right = 12  # cannot compare NamedArray and int

    err = "cannot compare <class 'ptools.namedarray.NamedArray'> and <class 'int'>"
    with pytest.raises(ValueError, match=err):
        assert left == right


def test_getitem():
    source = np.array((1, 2, 3, 4, 5))
    left = NamedArray("foo", "bar", source)
    assert_array_almost_equal(left.values, source)
    assert_array_almost_equal(left.values[1:4], source[1:4])
    assert_array_almost_equal(left[1:4].values, source[1:4])


def test_assignment():
    source = np.array((1, 2, 3, 4, 5))
    left = NamedArray("foo", "bar", source)
    right = np.array((6, 7, 8, 9, 10))
    left.values = right
    assert_array_almost_equal(left.values, right)

    left[0] = 42
    expected = np.array((42, 7, 8, 9, 10))
    assert_array_almost_equal(left.values, expected)


def test_copy():
    source = np.array((1, 2, 3, 4, 5))
    left = NamedArray("foo", "bar", source)
    right = left.copy()
    assert left == right
    assert left is not right

    # Check that the copy is a deep copy
    right.values[0] = 42
    assert_array_not_almost_equal(left.values, right.values)

    right.plural = "baz"
    assert left.plural != right.plural

    right.singular = "bat"
    assert left.singular != right.singular


def test_string_arrays_are_object():
    """Tests that string arrays are of dtype object.

    If not so, then the array will be truncated to the length of the
    shortest string in the array.
    """
    original = NamedArray("foo", "bar", np.array(("apple", "banana", "cherry")))
    original[2] = "oranges"
    assert np.all(original.values == np.array(("apple", "banana", "oranges")))
