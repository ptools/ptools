"""test_spatial - Tests for `ptools.spatial` module."""

# Unit-test specific imports.
import unittest
from .testing.moreassert import assert_array_almost_equal

# Scientific libraries.
import numpy as np

# Ptools imports.
from ptools import spatial
from ptools import linalg
from ptools.linalg import transform


# Ignores R0201: Method could be a function (no-self-use)
# pylint: disable=R0201
class TestSpatialObjectVector(unittest.TestCase):
    """Test spatial.SpatialObject methods on a vector."""

    def test_empty_constructor(self):
        obj = spatial.ObjectWithCoordinates()
        assert_array_almost_equal(obj.coords, (0, 0, 0))

    def test_constructor(self):
        obj = spatial.ObjectWithCoordinates((1, 1, 1))
        assert_array_almost_equal(obj.coords, (1, 1, 1))

    def test_translate_scalar(self):
        obj = spatial.TranslatableObject((0, 0, 0))
        obj.translate(1)
        assert_array_almost_equal(obj.coords, (1, 1, 1))

    def test_translate_vector(self):
        obj = spatial.TranslatableObject((0, 0, 0))
        obj.translate((1, 2, 3))
        assert_array_almost_equal(obj.coords, (1, 2, 3))

    def test_center_to_origin(self):
        obj = spatial.TranslatableObject((1, 1, 1))
        obj.center_to_origin()
        assert_array_almost_equal(obj.centroid(), (0, 0, 0))

    def test_center_custom_origin(self):
        obj = spatial.TranslatableObject((1, 1, 1))
        obj.center_to_origin(origin=(2, 2, 2))
        assert_array_almost_equal(obj.centroid(), (2, 2, 2))

    # pylint guesses type wrong
    # E1137: 'target.coords' does not support item assignment
    # pylint: disable=E1137
    def test_copy(self):
        source = spatial.ObjectWithCoordinates((6, 9, 12))
        target = source.copy()
        assert_array_almost_equal(target.coords, source.coords)
        target.coords[0] = 0
        assert_array_almost_equal(target.coords, [0, 9, 12])
        assert_array_almost_equal(source.coords, [6, 9, 12])

    def test_distance_to_axis(self):
        obj = spatial.ObjectWithCoordinates((0, 0, 0))
        axis = np.array((1, 0, 0))
        self.assertAlmostEqual(obj.distance_to_axis(axis), 0.0)

        obj.coords = (0, 1, 0)
        self.assertAlmostEqual(obj.distance_to_axis(axis), 1.0)

        obj.coords = (0, 1, 1)
        self.assertAlmostEqual(obj.distance_to_axis(axis), 2.0**0.5)


class TestSpatialObjectArray(unittest.TestCase):
    """Test spatial.SpatialObject methods on an array of vector."""

    def setUp(self):
        array = np.array(((0, 0, 0), (1, 1, 1)), dtype=float)
        self.obj = spatial.ObjectWithCoordinates(array)

    def test_constructor(self):
        self.assertEqual(self.obj.coords.shape, (2, 3))
        assert_array_almost_equal(self.obj.coords[1], (1, 1, 1))


class TestCoord3D(unittest.TestCase):
    def test_initialization_default(self):
        c = spatial.coord3d()
        self.assertEqual(c.shape, (3,))
        assert_array_almost_equal(c, (0, 0, 0))

    def test_initialization_scalar(self):
        c = spatial.coord3d(1)
        self.assertEqual(c.shape, (3,))
        assert_array_almost_equal(c, (1, 1, 1))

    def test_initialization_vector(self):
        c = spatial.coord3d([1, 2, 3])
        self.assertEqual(c.shape, (3,))
        assert_array_almost_equal(c, (1, 2, 3))

    def test_initialization_with_wrong_number_of_arguments(self):
        err = "Coordinates must be initialized either with 1 or 3 arguments"
        with self.assertRaisesRegex(ValueError, err):
            spatial.coord3d(12, 1)

        with self.assertRaisesRegex(ValueError, err):
            spatial.coord3d(12, 1, 2, 3)

        with self.assertRaisesRegex(ValueError, err):
            spatial.coord3d((1, 2), 3)

    def test_initialization_bad_dimensions(self):
        err = "3D coordinate array should be N x 3"
        with self.assertRaisesRegex(ValueError, err):
            spatial.coord3d([2, 2])

    def test_initialization_array(self):
        array = ((0, 0, 0), (1, 1, 1))
        c = spatial.coord3d(array)
        self.assertEqual(c.shape, (2, 3))
        assert_array_almost_equal(c[0], (0, 0, 0))
        assert_array_almost_equal(c[1], (1, 1, 1))

    def test_initialization_array_bad_dimensions(self):
        array = ((0, 0, 0, 0), (1, 1, 1, 1))
        err = "3D coordinate array should be N x 3"
        with self.assertRaisesRegex(ValueError, err):
            spatial.coord3d(array)

    def test_initialization_array_fails_more_than_2_dimensions(self):
        """Test raises the appropriate exeception when array has more that
        2 dimensions."""
        array = (
            ((0, 0, 0, 0), (1, 1, 1, 1)),
            ((0, 0, 0, 0), (1, 1, 1, 1)),
            ((0, 0, 0, 0), (1, 1, 1, 1)),
        )
        err = "3D coordinate array should have at most 2 dimensions"
        with self.assertRaisesRegex(ValueError, err):
            spatial.coord3d(array)


class TestSpatialObjectTransformations(unittest.TestCase):
    class FullTransformableObject(spatial.TranslatableObject, spatial.RotatableObject):
        pass

    def setUp(self):
        array = np.array(((0, 0, 0), (1, 1, 1)), dtype=float)
        self.obj = self.FullTransformableObject(array)

    def test_rotate(self):
        ref = [[0.0, 0.0, 0.0], [1.0, 0.8111595511436462, 1.1584558486938477]]
        matrix = linalg.rotation_matrix([10, 0, 0])
        self.obj.rotate(matrix)
        assert_array_almost_equal(self.obj.coords, ref)

    def test_rotate_with_vector(self):
        ref = [[0.0, 0.0, 0.0], [1.0, 0.8111595511436462, 1.1584558486938477]]
        transform.rotate(self.obj.coords, (10, 0, 0))
        assert_array_almost_equal(self.obj.coords, ref)

    def test_rotate_with_bad_input(self):
        err = r"Dimensions error: expected 3 x 1 or 4 x 4 \(got \(4,\)\)"
        with self.assertRaisesRegex(ValueError, err):
            self.obj.rotate((1, 0, 0, 0))

    def test_translate(self):
        self.obj.translate((1, 0, 0))
        assert_array_almost_equal(self.obj.coords, [[1, 0, 0], [2, 1, 1]])

    def test_translate_with_4_by_4(self):
        m = linalg.transformation_matrix(translation=(1, 0, 0))
        self.assertEqual(m.shape, (4, 4))
        self.obj.translate(m)
        assert_array_almost_equal(self.obj.coords, [[1, 0, 0], [2, 1, 1]])

    def test_translate_with_impossible_shape(self):
        err = r"Dimensions error: expected 3 x 1 or 4 x 4 \(got \(4,\)\)"
        with self.assertRaisesRegex(ValueError, err):
            self.obj.translate((1, 0, 0, 0))
