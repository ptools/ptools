""""test_atompropertycontainerlection - Tests for ``NamedArrayContainer``."""

from ptools.namedarray import NamedArray, NamedArrayContainer

import numpy as np
import pytest

from .testing import assert_array_almost_equal, assert_array_not_almost_equal


def generate_arrays() -> list[NamedArray]:
    return [
        NamedArray("one", "ones", np.ones(5)),
        NamedArray("two", "twos", np.ones(5) + 1),
        NamedArray("three", "threes", np.ones(5) + 3),
    ]

def test_initialization_empty():
    container = NamedArrayContainer()
    assert container.number_of_properties() == 0


def test_initialization_with_properties():
    property_list = generate_arrays()
    container = NamedArrayContainer(property_list)
    assert container.number_of_properties() == len(property_list)
    for prop in property_list:
        assert prop in container


def test_initialization_fails():
    properties = generate_arrays()
    properties[0].values = np.ones(properties[0].size() + 1)
    with pytest.raises(ValueError):
        NamedArrayContainer(properties)


def test_add_array():
    class expected:
        singular = "index"
        plural = "indices"
        values = (1, 2, 3)

    # Initialization.
    container = NamedArrayContainer()
    assert container.number_of_properties() == 0

    # Add a property and checks everything went ok.
    container.add_array(expected.singular, expected.plural, expected.values)
    assert container.number_of_properties() == 1
    assert container.number_of_elements() == 3
    assert expected.plural in container

    # Retrieves the property and checks its attributes.
    prop = container.get(expected.plural)
    assert prop.singular == expected.singular
    assert prop.plural == expected.plural
    assert_array_almost_equal(prop.values, expected.values)


def test_register():
    expected = NamedArray("index", "indices", (1, 2, 3))

    container = NamedArrayContainer()
    container.register(expected)

    actual = container.get(expected.plural)
    assert actual == expected


def test_register_checks_for_size():
    ones = NamedArray("one", "ones", np.ones(5))
    twos = NamedArray("two", "twos", np.ones(6) + 1)
    container = NamedArrayContainer()
    container.register(ones)
    with pytest.raises(ValueError):
        container.register(twos)


def test_register_copies_values():
    item = NamedArray("index", "indices", (1, 2, 3))
    container = NamedArrayContainer()
    container.register(item)

    assert_array_almost_equal(item.values, container.get(item.plural).values)

    item.values = 1
    assert_array_not_almost_equal(item.values, container.get(item.plural).values)
    assert_array_almost_equal((1, 2, 3), container.get(item.plural).values)


def test_getitem():
    container = NamedArrayContainer(generate_arrays())
    subset = container[1: 3]
    assert container.number_of_properties() == subset.number_of_properties()
    assert container.number_of_elements() == 5
    assert subset.number_of_elements() == 2

    subset = container[1]
    assert container.number_of_properties() == subset.number_of_properties()
    assert container.number_of_elements() == 5
    assert subset.number_of_elements() == 1


def test_equality():
    # Checks equality of containers containing identical arrays.
    lhs = NamedArrayContainer(generate_arrays())
    rhs = NamedArrayContainer(generate_arrays())
    assert lhs == rhs

    # Checks inequality of containers where one array differs from one container
    # to the other.
    r_ones = rhs.get("ones")
    r_ones.values[:] = 13
    assert lhs != rhs


def test_equality_fails():
    lhs = NamedArrayContainer(generate_arrays())
    rhs = 2
    with pytest.raises(TypeError):
        lhs == rhs


def test_copy():
    container = NamedArrayContainer(generate_arrays())
    copy = container.copy()
    assert container == copy

    # Checks that the copy is independent from the original.
    copy.get("ones").values[:] = 13
    assert container != copy


def test_add():
    lhs = NamedArrayContainer(generate_arrays())
    rhs = NamedArrayContainer(generate_arrays())
    result = lhs + rhs

    assert result.number_of_properties() == lhs.number_of_properties()
    assert result.number_of_elements() == lhs.number_of_elements() * 2
    for prop in lhs:
        assert_array_almost_equal(
            result.get(prop.plural).values,
            np.concatenate((lhs.get(prop.plural).values, rhs.get(prop.plural).values)),
        )


def test_set_values():
    container = NamedArrayContainer(generate_arrays())
    assert container.get("ones") == [1, 1, 1, 1, 1]

    container.set("ones", [2, 2, 2, 2, 2])
    assert container.get("ones") == [2, 2, 2, 2, 2]
