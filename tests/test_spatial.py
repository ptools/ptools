
"""test_spatial - Tests the `pyptools.spatial` module."""

import unittest

import pytest

from pyptools.spatial import SpatialObject, coord3d

from . import assert_array_almost_equal


class TestSpatialObject(unittest.TestCase):

    def test_empty_constructor(self):
        o = SpatialObject()
        assert_array_almost_equal(o.coords, (0, 0, 0))

    # def test_translate_vector(self):
    #     o = SpatialObject()



def test_coord3d_no_args():
    c = coord3d()
    assert c.shape == (3, )
    assert c[0] == c[1] == c[2] == 0


def test_coords3d_scalar():
    c = coord3d(1)
    assert c.shape == (3, )
    assert c[0] == c[1] == c[2] == 1

def test_coords3d_array():
    c = coord3d([2, 2, 2])
    assert c.shape == (3, )
    assert c[0] == c[1] == c[2] == 2

def test_coords3d_fails():
    with pytest.raises(ValueError) as excinfo:
        c = coord3d([2, 2])
    err = '3-d coordinates should be a scalar or 1 x 3 shaped-array'
    assert err in str(excinfo.value)
