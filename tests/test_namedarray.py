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

