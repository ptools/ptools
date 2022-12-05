""""test_atomproperty - Tests for ``AtomProperty``."""

from ptools.atom_property import AtomProperty

import numpy as np
import pytest

from .testing import assert_array_almost_equal, assert_array_not_almost_equal


def test_initialization_converts_to_numpy():
    prop = AtomProperty("foo", "bar", (1, 2, 3))
    assert isinstance(prop.values, np.ndarray)


def test_assignment_converts_to_numpy():
    prop = AtomProperty("foo", "bar", np.ones(5))
    prop.values = (1, 2, 3)
    assert isinstance(prop.values, np.ndarray)

    prop.values = 1
    assert isinstance(prop.values, np.ndarray)


def test_assignment_does_actual_copies():
    prop = AtomProperty("foo", "bar", np.ones(5))
    expected = np.array((1, 2, 3))

    prop.values = expected
    assert_array_almost_equal(prop.values, expected)

    expected[0] = 42
    expected[1] = 17
    expected[2] = 256

    assert_array_not_almost_equal(prop.values, expected)


def test_equality():
    left = AtomProperty("foo", "bar", np.ones(5))
    right = AtomProperty("foo", "bar", np.ones(5))
    assert left == right

    right = AtomProperty("foo", "bar", np.zeros(5))
    assert left != right


def test_equality_float():
    left = AtomProperty("foo", "bar", np.ones(5))
    right = AtomProperty("foo", "bar", np.ones(5))
    assert left == right

    right = AtomProperty("foo", "bar", np.ones(5) + 1e-6)
    assert left == right

    right = AtomProperty("foo", "bar", np.ones(5) + 1e-4)
    assert left != right

    # Check array comparison function works as expected
    left._array_comparison_func = lambda x, y: True
    assert left == right


def test_equality_type_error():
    left = AtomProperty("foo", "bar", np.ones(5))
    right = 12  # cannot compare AtomProperty and int

    err = "cannot compare <class 'ptools.atom_property.AtomProperty'> and <class 'int'>"
    with pytest.raises(ValueError, match=err):
        assert left == right


def test_getitem():
    source = np.array((1, 2, 3, 4, 5))
    left = AtomProperty("foo", "bar", source)
    assert_array_almost_equal(left.values, source)
    assert_array_almost_equal(left.values[1:4], source[1:4])
    assert_array_almost_equal(left[1:4].values, source[1:4])

