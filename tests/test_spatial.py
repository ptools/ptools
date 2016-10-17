
"""test_spatial - Tests the `pyptools.spatial` module."""

import unittest

import pytest

from pyptools.spatial import SpatialObject, coord3d

from . import assert_array_almost_equal


class TestSpatialObjectVector(unittest.TestCase):
    """Test SpatialObject methods on a vector."""

    def test_empty_constructor(self):
        o = SpatialObject()
        assert_array_almost_equal(o.coords, (0, 0, 0))

    def test_constructor(self):
        o = SpatialObject((1, 1, 1))
        assert_array_almost_equal(o.coords, (1, 1, 1))

    def test_translate_scalar(self):
        o = SpatialObject((0, 0, 0))
        o.translate(1)
        assert_array_almost_equal(o.coords, (1, 1, 1))

    def test_translate_vector(self):
        o = SpatialObject((0, 0, 0))
        o.translate((1, 2, 3))
        assert_array_almost_equal(o.coords, (1, 2, 3))


class TestSpatialObjectArray(unittest.TestCase):
    """Test SpatialObject methods on an array of vector."""

    def setUp(self):
        import numpy
        array = numpy.array(((0, 0, 0), (1, 1, 1)), dtype=float)
        self.obj = SpatialObject(array)

    def test_constructor(self):
        self.assertEqual(self.obj.coords.shape, (2, 3))
        assert_array_almost_equal(self.obj.coords[1], (1, 1, 1))


def test_coord3d_no_args():
    c = coord3d()
    assert c.shape == (3, )
    assert_array_almost_equal(c, (0, 0, 0))


def test_coord3d_scalar():
    c = coord3d(1)
    assert c.shape == (3, )
    assert_array_almost_equal(c, (1, 1, 1))


def test_coord3d_vector():
    c = coord3d([2, 2, 2])
    assert c.shape == (3, )
    assert_array_almost_equal(c, (2, 2, 2))


def test_coord3d_vector_fails():
    with pytest.raises(ValueError) as excinfo:
        coord3d([2, 2])
    err = '3-d coordinates should be a scalar or 1 x 3 shaped-array'
    assert err in str(excinfo.value)


def test_coord3d_array():
    array = ((0, 0, 0), (1, 1, 1))
    c = coord3d(array)
    assert c.shape == (2, 3)
    assert_array_almost_equal(c[0], (0, 0, 0))
    assert_array_almost_equal(c[1], (1, 1, 1))


def test_coord3d_array_fails():
    array = ((0, 0, 0, 0), (1, 1, 1, 1))
    with pytest.raises(ValueError) as excinfo:
        coord3d(array)
    err = '3-d coordinate array should be N x 3'
    assert err in str(excinfo.value)


def test_coord3d_array_fails2():
    """Test raises the appropriate exeception when array has more that
    2 dimensions."""
    array = (((0, 0, 0, 0), (1, 1, 1, 1)),
             ((0, 0, 0, 0), (1, 1, 1, 1)),
             ((0, 0, 0, 0), (1, 1, 1, 1)))
    with pytest.raises(ValueError) as excinfo:
        coord3d(array)
    err = '3-d coordinate array should have at most 2 dimensions'
    assert err in str(excinfo.value)
