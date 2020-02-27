
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


class TestRotation(unittest.TestCase):

    def test_rotation_X(self):
        # Rotation by 90°.
        ref = [[  1.00,    0.00,    0.00,    0.00],
               [  0.00,    0.00,   -1.00,    0.00],
               [  0.00,    1.00,    0.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation_by(alpha=90)
        assert_array_almost_equal(m, ref)

        # Rotation by 180°.
        ref = [[  1.00,    0.00,    0.00,    0.00],
               [  0.00,   -1.00,   -0.00,    0.00],
               [  0.00,    0.00,   -1.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation_by(alpha=180)
        assert_array_almost_equal(m, ref)

        # Rotation by 10°.
        ref = [[ 1.0000000,    0.000000,    0.000000,    0.000000],
               [ 0.0000000,    0.984808,   -0.173648,    0.000000],
               [ 0.0000000,    0.173648,    0.984808,    0.000000],
               [ 0.0000000,    0.000000,    0.000000,    1.000000]]
        m = spatial.rotation_by(alpha=10)
        assert_array_almost_equal(m, ref)

    def test_rotation_Y(self):
        # Rotation by 90°.
        ref = [[  0.00,    0.00,    1.00,    0.00],
               [  0.00,    1.00,    0.00,    0.00],
               [ -1.00,    0.00,    0.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation_by(beta=90)
        assert_array_almost_equal(m, ref)

        # Rotation by 180°.
        ref = [[ -1.00,    0.00,    0.00,    0.00],
               [  0.00,    1.00,    0.00,    0.00],
               [ -0.00,    0.00,   -1.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation_by(beta=180)
        assert_array_almost_equal(m, ref)

        # Rotation by 10°.
        ref = [[ 0.9848078,    0.000000,    0.173648,    0.000000],
               [ 0.0000000,    1.000000,    0.000000,    0.000000],
               [-0.1736482,    0.000000,    0.984808,    0.000000],
               [ 0.0000000,    0.000000,    0.000000,    1.000000]]
        m = spatial.rotation_by(beta=10)
        assert_array_almost_equal(m, ref)

    def test_rotation_Z(self):
        # Rotation by 90°.
        ref = [[  0.00,   -1.00,    0.00,    0.00],
               [  1.00,    0.00,    0.00,    0.00],
               [  0.00,    0.00,    1.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation_by(gamma=90)
        assert_array_almost_equal(m, ref)

        # Rotation by 180°.
        ref = [[ -1.00,   -0.00,    0.00,    0.00],
               [  0.00,   -1.00,    0.00,    0.00],
               [  0.00,    0.00,    1.00,    0.00],
               [  0.00,    0.00,    0.00,    1.00]]
        m = spatial.rotation_by(gamma=180)
        assert_array_almost_equal(m, ref)

        # Rotation by 10°.
        ref = [[ 0.9848078,   -0.173648,    0.000000,    0.000000],
               [ 0.1736482,    0.984808,    0.000000,    0.000000],
               [ 0.0000000,    0.000000,    1.000000,    0.000000],
               [ 0.0000000,    0.000000,    0.000000,    1.000000]]
        m = spatial.rotation_by(gamma=10)
        assert_array_almost_equal(m, ref)

    def test_rotation_XYZ(self):
        """Test rotation in X, Y and Z is equivalent to successive rotation
        in X, then Y, then Z."""
        m = spatial.rotation_by(alpha=10, beta=10, gamma=10)
        ref = spatial.rotation_by(alpha=10).dot(spatial.rotation_by(beta=10)).dot(spatial.rotation_by(gamma=10))
        assert_array_almost_equal(m, ref)

    def test_rotate_by(self):
        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = numpy.array(list([i + 1., i + 11., i + 21.] for i in range(10)))

        # Rotate by 12° along X-axis.
        spatial.rotate_by(coords, alpha=12)

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

    def test_rotate_spatial_object(self):
        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = numpy.array(list([i + 1., i + 11., i + 21.] for i in range(10)))

        # Initialize SpatialObject.
        o = spatial.SpatialObject(coords)

        # Rotate by 12° along X-axis.
        o.rotate_by(alpha=12)

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

    def test_rotate_abstract_euler(self):
        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = numpy.array(list([i + 1., i + 11., i + 21.] for i in range(10)))

        # Initialize SpatialObject.
        o = spatial.SpatialObject(coords)

        # Attract Euler rotation.
        o.attract_euler_rotate(10, 12, 14)

        # Coordinates calculated with ptools version c9f7fee::
        #
        #     >>> from ptools import *
        #     >>> rb = Rigidbody('test_10atoms.pdb')
        #     >>> rb.AttractEulerRotate(10, 12, 14)
        #     >>> for i in xrange(10):
        #     ...     c = r.getCoords(i)
        #     ...     print c.x, c.y, c.z
        #
        ref = numpy.array([
            [1.92178620509, 0.634022491719, 23.6411664954],
            [1.10926523817, 1.12485261069, 25.0899230217],
            [0.296744271247, 1.61568272966, 26.5386795481],
            [-0.515776695676, 2.10651284863, 27.9874360744],
            [-1.3282976626, 2.5973429676, 29.4361926007],
            [-2.14081862952, 3.08817308657, 30.8849491271],
            [-2.95333959645, 3.57900320554, 32.3337056534],
            [-3.76586056337, 4.06983332451, 33.7824621798],
            [-4.57838153029, 4.56066344348, 35.2312187061],
            [-5.39090249721, 5.05149356245, 36.6799752325],
        ])

        assert_array_almost_equal(o.coords, ref)


class TestTransformation(unittest.TestCase):
    """Test transformations thanks to 4 x 4 matrices."""

    def test_transform_rotation(self):
        """Test transformation with a rotation matrix.

        In this test, the last row and columns are identity.
        """
        # Matrix calculated with VMD that correspond to a rotation
        # in the X-axis by 10°.
        m = numpy.array([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.984807753012208, -0.17364817766693033, 0.0],
            [0.0, 0.17364817766693033, 0.984807753012208, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])

        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = numpy.array(list([i + 1., i + 11., i + 21.] for i in range(10)))

        ref = coords.copy()
        spatial.rotate_by(ref, alpha=10)

        spatial.transform(coords, m)
        assert_array_almost_equal(coords, ref)

    def test_transform(self):
        """Test transformation with a random matrix."""
        # Input transformation matrix (just a random matrix).
        m = numpy.array([[ 0.72898284,  0.43879558,  0.31011219,  0.32696626],
                         [ 0.74291592,  0.00244003,  0.20881428,  0.91014385],
                         [ 0.71045798,  0.85828462,  0.71570462,  0.96008097],
                         [ 0.43315630,  0.30949433,  0.93486660,  0.48177513]])

        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = numpy.array(list([i + 1., i + 11., i + 21.] for i in range(10)))

        # Reference calculated with VMD
        # with input topology = tests/data/test_10atoms.pdb.
        ref = numpy.array([
            [0.5175050497055054, 0.253219336271286, 1.0914303064346313],
            [0.5412969589233398, 0.2738751769065857, 1.1091270446777344],
            [0.5621657967567444, 0.291993111371994, 1.1246496438980103],
            [0.5806189775466919, 0.3080138564109802, 1.1383754014968872],
            [0.5970528721809387, 0.3222815990447998, 1.1505991220474243],
            [0.611781895160675, 0.33506909012794495, 1.1615549325942993],
            [0.6250582337379456, 0.3465954065322876, 1.171429991722107],
            [0.6370866894721985, 0.35703831911087036, 1.1803770065307617],
            [0.6480352282524109, 0.3665436804294586, 1.1885206699371338],
            [0.6580431461334229, 0.37523239850997925, 1.1959648132324219]
        ])

        spatial.transform(coords, m)
        assert_array_almost_equal(coords, ref)

    def test_transform_spatial_object(self):
        """Test spatial.SpatialObject.transform."""
        # Input transformation matrix (just a random matrix).
        m = numpy.array([[ 0.72898284,  0.43879558,  0.31011219,  0.32696626],
                         [ 0.74291592,  0.00244003,  0.20881428,  0.91014385],
                         [ 0.71045798,  0.85828462,  0.71570462,  0.96008097],
                         [ 0.43315630,  0.30949433,  0.93486660,  0.48177513]])

        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = numpy.array(list([i + 1., i + 11., i + 21.] for i in range(10)))

        # Initialize SpatialObject.
        o = spatial.SpatialObject(coords)

        # Rotate by 12° along X-axis.
        o.transform(m)

        # Reference calculated with VMD
        # with input topology = tests/data/test_10atoms.pdb.
        ref = numpy.array([
            [0.5175050497055054, 0.253219336271286, 1.0914303064346313],
            [0.5412969589233398, 0.2738751769065857, 1.1091270446777344],
            [0.5621657967567444, 0.291993111371994, 1.1246496438980103],
            [0.5806189775466919, 0.3080138564109802, 1.1383754014968872],
            [0.5970528721809387, 0.3222815990447998, 1.1505991220474243],
            [0.611781895160675, 0.33506909012794495, 1.1615549325942993],
            [0.6250582337379456, 0.3465954065322876, 1.171429991722107],
            [0.6370866894721985, 0.35703831911087036, 1.1803770065307617],
            [0.6480352282524109, 0.3665436804294586, 1.1885206699371338],
            [0.6580431461334229, 0.37523239850997925, 1.1959648132324219]
        ])

        assert_array_almost_equal(o.coords, ref)


class TestCoord3D(unittest.TestCase):

    def test_empty_constructor(self):
        c = spatial.coord3d()
        self.assertEqual(c.shape, (3, ))
        assert_array_almost_equal(c, (0, 0, 0))

    def test_constructor_scalar(self):
        c = spatial.coord3d(1)
        self.assertEqual(c.shape, (3, ))
        assert_array_almost_equal(c, (1, 1, 1))

    def test_constructor_vector(self):
        c = spatial.coord3d([1, 2, 3])
        self.assertEqual(c.shape, (3, ))
        assert_array_almost_equal(c, (1, 2, 3))

    def test_constructor_vector_fails(self):
        err = '3-d coordinates should be a scalar or 1 x 3 shaped-array'
        with self.assertRaisesRegex(ValueError, err):
            spatial.coord3d([2, 2])

    def test_constructor_array(self):
        array = ((0, 0, 0), (1, 1, 1))
        c = spatial.coord3d(array)
        self.assertEqual(c.shape, (2, 3))
        assert_array_almost_equal(c[0], (0, 0, 0))
        assert_array_almost_equal(c[1], (1, 1, 1))

    def test_constructor_array_fails(self):
        array = ((0, 0, 0, 0), (1, 1, 1, 1))
        err = '3-d coordinate array should be N x 3'
        with self.assertRaisesRegex(ValueError, err):
            spatial.coord3d(array)

    def test_constructor_array_fails2(self):
        """Test raises the appropriate exeception when array has more that
        2 dimensions."""
        array = (((0, 0, 0, 0), (1, 1, 1, 1)),
                 ((0, 0, 0, 0), (1, 1, 1, 1)),
                 ((0, 0, 0, 0), (1, 1, 1, 1)))
        err = '3-d coordinate array should have at most 2 dimensions'
        with self.assertRaisesRegex(ValueError, err):
            spatial.coord3d(array)


class TestSpatialObjectTransformations(unittest.TestCase):
    def setUp(self):
        array = numpy.array(((0, 0, 0), (1, 1, 1)), dtype=float)
        self.obj = spatial.SpatialObject(array)

    def test_rotate(self):
        ref = [[0.0, 0.0, 0.0],
               [1.0, 0.8111595511436462, 1.1584558486938477]]
        matrix = spatial.rotation_by(alpha=10)
        self.obj.rotate(matrix)
        assert_array_almost_equal(self.obj.coords, ref)
