""""test_atompropertycollection - Tests for ``AtomPropertyCollection``."""

from ptools.atom_property import AtomProperty, AtomPropertyCollection


from .testing import assert_array_almost_equal, assert_array_not_almost_equal
import pytest


def test_initialization_empty():
    col = AtomPropertyCollection()
    assert len(col) == 0


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


def test_register_copies_values():
    item = AtomProperty("index", "indices", (1, 2, 3))
    col = AtomPropertyCollection()
    col.register(item)

    assert_array_almost_equal(item.values, col[item.plural].values)

    item.values = 1
    assert_array_not_almost_equal(item.values, col[item.plural].values)
    assert_array_almost_equal((1, 2, 3), col[item.plural].values)
