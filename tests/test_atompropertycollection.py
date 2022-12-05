""""test_atompropertycollection - Tests for ``AtomPropertyCollection``."""

from ptools.atom_property import AtomProperty, AtomPropertyCollection

import numpy as np
import pytest

from .testing import assert_array_almost_equal, assert_array_not_almost_equal


def generate_properties():
    return [
        AtomProperty("one", "ones", np.ones(5)),
        AtomProperty("two", "twos", np.ones(5) + 1),
        AtomProperty("three", "threes", np.ones(5) + 3),
    ]

def test_initialization_empty():
    col = AtomPropertyCollection()
    assert len(col) == 0


def test_initialization_with_properties():
    property_list = generate_properties()
    col = AtomPropertyCollection(property_list)
    assert len(col) == len(property_list)
    for prop in property_list:
        assert prop in col


def test_initialization_fails():
    properties = generate_properties()
    properties[0].values = np.ones(properties[0].size() + 1)
    with pytest.raises(ValueError):
        AtomPropertyCollection(properties)


def test_add_property():
    singular = "index"
    plural = "indices"
    values = (1, 2, 3)

    col = AtomPropertyCollection()
    col.add_property(singular, plural, values)

    assert len(col) == 1
    assert plural in col
    prop = col[plural]
    assert prop.singular == singular
    assert prop.plural == plural
    assert_array_almost_equal(values, prop.values)


def test_register():
    expected = AtomProperty("index", "indices", (1, 2, 3))

    col = AtomPropertyCollection()
    col.register(expected)

    actual = col[expected.plural]
    assert actual == expected


def test_register_checks_for_size():
    ones = AtomProperty("one", "ones", np.ones(5))
    twos = AtomProperty("two", "twos", np.ones(6) + 1)
    col = AtomPropertyCollection()
    col.register(ones)
    with pytest.raises(ValueError):
        col.register(twos)


def test_register_copies_values():
    item = AtomProperty("index", "indices", (1, 2, 3))
    col = AtomPropertyCollection()
    col.register(item)

    assert_array_almost_equal(item.values, col[item.plural].values)

    item.values = 1
    assert_array_not_almost_equal(item.values, col[item.plural].values)
    assert_array_almost_equal((1, 2, 3), col[item.plural].values)


def test_slice():
    collection = AtomPropertyCollection(generate_properties())
    subset = collection.slice(slice(1, 3))
    assert collection.number_of_properties() == subset.number_of_properties()
    assert collection.number_of_elements() == 5
    assert subset.number_of_elements() == 2

    subset = collection.slice(1)
    assert collection.number_of_properties() == subset.number_of_properties()
    assert collection.number_of_elements() == 5
    assert subset.number_of_elements() == 1

