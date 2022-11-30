"""test_spatial - Tests for `ptools.spatial` module."""

# Unit-test specific imports.
import unittest
from .testing.moreassert import assert_array_almost_equal

# Scientific libraries.
import numpy as np

# Ptools imports.
from ptools import linalg, measure, spatial
from ptools.linalg import transform


# Ignores R0201: Method could be a function (no-self-use)
# pylint: disable=R0201
class TestSpatialObjectVector(unittest.TestCase):
    """Test spatial.SpatialObject methods on a vector."""

    def test_empty_constructor(self):
        obj = spatial.ObjectWithCoordinates()
        assert_array_almost_equal(obj.coordinates, (0, 0, 0))

    def test_constructor(self):
        obj = spatial.ObjectWithCoordinates((1, 1, 1))
        assert_array_almost_equal(obj.coordinates, (1, 1, 1))

    def test_translate_scalar(self):
        obj = spatial.SupportsTranslation((0, 0, 0))
        obj.translate(1)
        assert_array_almost_equal(obj.coordinates, (1, 1, 1))

    def test_translate_vector(self):
        obj = spatial.SupportsTranslation((0, 0, 0))
        obj.translate((1, 2, 3))
        assert_array_almost_equal(obj.coordinates, (1, 2, 3))

    def test_center_to_origin(self):
        obj = spatial.SupportsTranslation((1, 1, 1))
        obj.center_to_origin()
        assert_array_almost_equal(measure.centroid(obj), (0, 0, 0))

    def test_center_custom_origin(self):
        obj = spatial.SupportsTranslation((1, 1, 1))
        obj.center_to_origin(origin=(2, 2, 2))
        assert_array_almost_equal(measure.centroid(obj), (2, 2, 2))

    # pylint guesses type wrong
    # E1137: 'target.coordinates' does not support item assignment
    # pylint: disable=E1137
    def test_copy(self):
        source = spatial.ObjectWithCoordinates((6, 9, 12))
        target = source.copy()
        assert_array_almost_equal(target.coordinates, source.coordinates)
        target.coordinates[0] = 0
        assert_array_almost_equal(target.coordinates, [0, 9, 12])
        assert_array_almost_equal(source.coordinates, [6, 9, 12])


class TestSpatialObjectArray(unittest.TestCase):
    """Test spatial.SpatialObject methods on an array of vector."""

    def setUp(self):
        array = np.array(((0, 0, 0), (1, 1, 1)), dtype=float)
        self.obj = spatial.ObjectWithCoordinates(array)

    def test_constructor(self):
        self.assertEqual(self.obj.coordinates.shape, (2, 3))
        assert_array_almost_equal(self.obj.coordinates[1], (1, 1, 1))


class TestSpatialObjectTransformations(unittest.TestCase):
    class FullTransformableObject(spatial.SupportsTranslation, spatial.SupportsRotation):
        pass

    def setUp(self):
        array = np.array(((0, 0, 0), (1, 1, 1)), dtype=float)
        self.obj = self.FullTransformableObject(array)

    def test_rotate(self):
        ref = [[0.0, 0.0, 0.0], [1.0, 0.8111595511436462, 1.1584558486938477]]
        matrix = linalg.rotation_matrix([10, 0, 0])
        self.obj.rotate(matrix)
        assert_array_almost_equal(self.obj.coordinates, ref)

    def test_rotate_with_vector(self):
        ref = [[0.0, 0.0, 0.0], [1.0, 0.8111595511436462, 1.1584558486938477]]
        transform.rotate(self.obj.coordinates, (10, 0, 0))
        assert_array_almost_equal(self.obj.coordinates, ref)

    def test_rotate_with_bad_input(self):
        err = r"Dimensions error: expected 3 x 1 or 4 x 4 \(got \(4,\)\)"
        with self.assertRaisesRegex(ValueError, err):
            self.obj.rotate((1, 0, 0, 0))

    def test_translate(self):
        self.obj.translate((1, 0, 0))
        assert_array_almost_equal(self.obj.coordinates, [[1, 0, 0], [2, 1, 1]])

    def test_translate_with_4_by_4(self):
        m = linalg.transformation_matrix(translation=(1, 0, 0))
        self.assertEqual(m.shape, (4, 4))
        self.obj.translate(m)
        assert_array_almost_equal(self.obj.coordinates, [[1, 0, 0], [2, 1, 1]])

    def test_translate_with_impossible_shape(self):
        err = r"Dimensions error: expected 3 x 1 or 4 x 4 \(got \(4,\)\)"
        with self.assertRaisesRegex(ValueError, err):
            self.obj.translate((1, 0, 0, 0))
