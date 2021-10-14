"""testing.moreassert - Defines more assert functions."""
import unittest

import numpy

from ptools.atom import BaseAtom

from .dummy import DUMMY_ATOM_ATTRS


assert_array_almost_equal = numpy.testing.assert_array_almost_equal
assert_array_equal = numpy.testing.assert_array_equal


def assert_array_not_equal(array1, array2):
    """Asserts numpy arrays are not equal."""
    try:
        assert_array_equal(array1, array2)
    except AssertionError:
        return
    raise AssertionError("arrays are equal")


def assert_array_not_almost_equal(array1, array2):
    """Asserts numpy arrays are not equal."""
    try:
        assert_array_almost_equal(array1, array2)
    except AssertionError:
        return
    raise AssertionError("arrays are almost equal")


def assert_dummy_atom_initialization_ok(self: unittest.TestCase, dummy: BaseAtom):
    """Asserts dummy atom initialization went fine."""
    for attr, value in DUMMY_ATOM_ATTRS.items():
        if attr != "coords":
            self.assertEqual(getattr(dummy, attr), value)
    assert_array_almost_equal(dummy.coords, DUMMY_ATOM_ATTRS["coords"])
