
"""test_spatial - Tests for `pyptools.spatial` module."""

import unittest

import numpy
import pytest


import pyptools.spatial as spatial


from .testing.moreassert import assert_array_almost_equal


class TestSpatialObjectVector(unittest.TestCase):
    """Test spatial.SpatialObject methods on a vector."""

    def test_empty_constructor(self):
        o = spatial.SpatialObject()
        assert_array_almost_equal(o.coords, (0, 0, 0))

    def test_constructor(self):
        o = spatial.SpatialObject((1, 1, 1))
        assert_array_almost_equal(o.coords, (1, 1, 1))

    def test_translate_scalar(self):
        o = spatial.SpatialObject((0, 0, 0))
        o.translate(1)
        assert_array_almost_equal(o.coords, (1, 1, 1))

    def test_translate_vector(self):
        o = spatial.SpatialObject((0, 0, 0))
        o.translate((1, 2, 3))
        assert_array_almost_equal(o.coords, (1, 2, 3))


class TestSpatialObjectArray(unittest.TestCase):
    """Test spatial.SpatialObject methods on an array of vector."""

    def setUp(self):
        array = numpy.array(((0, 0, 0), (1, 1, 1)), dtype=float)
        self.obj = spatial.SpatialObject(array)

    def test_constructor(self):
        self.assertEqual(self.obj.coords.shape, (2, 3))
        assert_array_almost_equal(self.obj.coords[1], (1, 1, 1))

    def test_rotate(self):
        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = numpy.array(list([i + 1., i + 11., i + 21.] for i in range(10)))

        # Initialize SpatialObject.
        o = spatial.SpatialObject(coords)

        # Rotate by 12° along X-axis.
        o.rotate(alpha=12)

        # Reference coordinates calculated with VMD::
        #     set m [transaxis x 12]
        #     set sel [atomselect 0 "all"]
        #     $sel move $m
        ref_coords = [[ 1.0000000,    6.393478,   22.828129],
                      [ 2.0000000,    7.163715,   24.014187],
                      [ 3.0000000,    7.933951,   25.200247],
                      [ 4.0000000,    8.704186,   26.386307],
                      [ 5.0000000,    9.474422,   27.572367],
                      [ 6.0000000,   10.244658,   28.758427],
                      [ 7.0000000,   11.014894,   29.944485],
                      [ 8.0000000,   11.785130,   31.130545],
                      [ 9.0000000,   12.555366,   32.316605],
                      [10.0000000,   13.325602,   33.502663]]

        assert_array_almost_equal(o.coords, ref_coords, decimal=5)


class TestRotation(unittest.TestCase):

    def test_rotation_X(self):
        # Rotation by 90°.
        ref = [[  1.00,    0.00,    0.00,    0.00],
               [  0.00,    0.00,   -1.00,    0.00],
               [  0.00,    1.00,    0.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation(alpha=90)
        assert_array_almost_equal(m, ref)

        # Rotation by 180°.
        ref = [[  1.00,    0.00,    0.00,    0.00],
               [  0.00,   -1.00,   -0.00,    0.00],
               [  0.00,    0.00,   -1.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation(alpha=180)
        assert_array_almost_equal(m, ref)

        # Rotation by 10°.
        ref = [[ 1.0000000,    0.000000,    0.000000,    0.000000],
               [ 0.0000000,    0.984808,   -0.173648,    0.000000],
               [ 0.0000000,    0.173648,    0.984808,    0.000000],
               [ 0.0000000,    0.000000,    0.000000,    1.000000]]
        m = spatial.rotation(alpha=10)
        assert_array_almost_equal(m, ref)

    def test_rotation_Y(self):
        # Rotation by 90°.
        ref = [[  0.00,    0.00,    1.00,    0.00],
               [  0.00,    1.00,    0.00,    0.00],
               [ -1.00,    0.00,    0.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation(beta=90)
        assert_array_almost_equal(m, ref)

        # Rotation by 180°.
        ref = [[ -1.00,    0.00,    0.00,    0.00],
               [  0.00,    1.00,    0.00,    0.00],
               [ -0.00,    0.00,   -1.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation(beta=180)
        assert_array_almost_equal(m, ref)

        # Rotation by 10°.
        ref = [[ 0.9848078,    0.000000,    0.173648,    0.000000],
               [ 0.0000000,    1.000000,    0.000000,    0.000000],
               [-0.1736482,    0.000000,    0.984808,    0.000000],
               [ 0.0000000,    0.000000,    0.000000,    1.000000]]
        m = spatial.rotation(beta=10)
        assert_array_almost_equal(m, ref)

    def test_rotation_Z(self):
        # Rotation by 90°.
        ref = [[  0.00,   -1.00,    0.00,    0.00],
               [  1.00,    0.00,    0.00,    0.00],
               [  0.00,    0.00,    1.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation(gamma=90)
        assert_array_almost_equal(m, ref)

        # Rotation by 180°.
        ref = [[ -1.00,   -0.00,    0.00,    0.00],
               [  0.00,   -1.00,    0.00,    0.00],
               [  0.00,    0.00,    1.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation(gamma=180)
        assert_array_almost_equal(m, ref)

        # Rotation by 10°.
        ref = [[ 0.9848078,   -0.173648,    0.000000,    0.000000],
               [ 0.1736482,    0.984808,    0.000000,    0.000000],
               [ 0.0000000,    0.000000,    1.000000,    0.000000],
               [ 0.0000000,    0.000000,    0.000000,    1.000000]]
        m = spatial.rotation(gamma=10)
        assert_array_almost_equal(m, ref)

    def test_rotate(self):
        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = numpy.array(list([i + 1., i + 11., i + 21.] for i in range(10)))

        # Rotate by 12° along X-axis.
        spatial.rotate(coords, alpha=12)

        # Reference coordinates calculated with VMD::
        #     set m [transaxis x 12]
        #     set sel [atomselect 0 "all"]
        #     $sel move $m
        ref_coords = [[ 1.0000000,    6.393478,   22.828129],
                      [ 2.0000000,    7.163715,   24.014187],
                      [ 3.0000000,    7.933951,   25.200247],
                      [ 4.0000000,    8.704186,   26.386307],
                      [ 5.0000000,    9.474422,   27.572367],
                      [ 6.0000000,   10.244658,   28.758427],
                      [ 7.0000000,   11.014894,   29.944485],
                      [ 8.0000000,   11.785130,   31.130545],
                      [ 9.0000000,   12.555366,   32.316605],
                      [10.0000000,   13.325602,   33.502663]]

        assert_array_almost_equal(coords, ref_coords, decimal=5)


def test_coord3d_no_args():
    c = spatial.coord3d()
    assert c.shape == (3, )
    assert_array_almost_equal(c, (0, 0, 0))


def test_coord3d_scalar():
    c = spatial.coord3d(1)
    assert c.shape == (3, )
    assert_array_almost_equal(c, (1, 1, 1))


def test_coord3d_vector():
    c = spatial.coord3d([2, 2, 2])
    assert c.shape == (3, )
    assert_array_almost_equal(c, (2, 2, 2))


def test_coord3d_vector_fails():
    with pytest.raises(ValueError) as excinfo:
        spatial.coord3d([2, 2])
    err = '3-d coordinates should be a scalar or 1 x 3 shaped-array'
    assert err in str(excinfo.value)


def test_coord3d_array():
    array = ((0, 0, 0), (1, 1, 1))
    c = spatial.coord3d(array)
    assert c.shape == (2, 3)
    assert_array_almost_equal(c[0], (0, 0, 0))
    assert_array_almost_equal(c[1], (1, 1, 1))


def test_coord3d_array_fails():
    array = ((0, 0, 0, 0), (1, 1, 1, 1))
    with pytest.raises(ValueError) as excinfo:
        spatial.coord3d(array)
    err = '3-d coordinate array should be N x 3'
    assert err in str(excinfo.value)


def test_coord3d_array_fails2():
    """Test raises the appropriate exeception when array has more that
    2 dimensions."""
    array = (((0, 0, 0, 0), (1, 1, 1, 1)),
             ((0, 0, 0, 0), (1, 1, 1, 1)),
             ((0, 0, 0, 0), (1, 1, 1, 1)))
    with pytest.raises(ValueError) as excinfo:
        spatial.coord3d(array)
    err = '3-d coordinate array should have at most 2 dimensions'
    assert err in str(excinfo.value)
