"""testing.moreassert - Defines more assert functions."""
import unittest

import numpy

from ptools.atom import BaseAtom


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



def dummy_atom():
    """Creates a dummy atoms used in several tests."""
    return BaseAtom(
        name="CA",
        resname="ALA",
        chain="A",
        index=42,
        resid=17,
        charge=2.0,
        coords=(1, 2, 3),
    )

def assert_dummy_atom_initialization_ok(self: unittest.TestCase, dummy: BaseAtom):
    self.assertEqual(dummy.name, "CA")
    self.assertEqual(dummy.resname, "ALA")
    self.assertEqual(dummy.chain, "A")
    self.assertEqual(dummy.index, 42)
    self.assertEqual(dummy.resid, 17)
    self.assertEqual(dummy.charge, 2.0)
    assert_array_almost_equal(dummy.coords, (1, 2, 3))
