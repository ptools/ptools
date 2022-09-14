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
        obj = spatial.SpatialObject()
        assert_array_almost_equal(obj.coords, (0, 0, 0))

    def test_constructor(self):
        obj = spatial.SpatialObject((1, 1, 1))
        assert_array_almost_equal(obj.coords, (1, 1, 1))

    def test_translate_scalar(self):
        obj = spatial.SpatialObject((0, 0, 0))
        obj.translate(1)
        assert_array_almost_equal(obj.coords, (1, 1, 1))

    def test_translate_vector(self):
        obj = spatial.SpatialObject((0, 0, 0))
        obj.translate((1, 2, 3))
        assert_array_almost_equal(obj.coords, (1, 2, 3))

    def test_center_to_origin(self):
        obj = spatial.SpatialObject((1, 1, 1))
        obj.center_to_origin()
        assert_array_almost_equal(obj.centroid(), (0, 0, 0))

    def test_center_custom_origin(self):
        obj = spatial.SpatialObject((1, 1, 1))
        obj.center_to_origin(origin=(2, 2, 2))
        assert_array_almost_equal(obj.centroid(), (2, 2, 2))

    # pylint guesses type wrong
    # E1137: 'target.coords' does not support item assignment
    # pylint: disable=E1137
    def test_copy(self):
        source = spatial.SpatialObject((6, 9, 12))
        target = source.copy()
        assert_array_almost_equal(target.coords, source.coords)
        target.coords[0] = 0
        assert_array_almost_equal(target.coords, [0, 9, 12])
        assert_array_almost_equal(source.coords, [6, 9, 12])

    def test_distance_to_axis(self):
        obj = spatial.SpatialObject((0, 0, 0))
        axis = np.array((1, 0, 0))
        self.assertAlmostEqual(obj.distance_to_axis(axis), 0.0)

        obj.coords = (0, 1, 0)
        self.assertAlmostEqual(obj.distance_to_axis(axis), 1.0)

        obj.coords = (0, 1, 1)
        self.assertAlmostEqual(obj.distance_to_axis(axis), 2.0 ** 0.5)


class TestSpatialObjectArray(unittest.TestCase):
    """Test spatial.SpatialObject methods on an array of vector."""

    def setUp(self):
        array = np.array(((0, 0, 0), (1, 1, 1)), dtype=float)
        self.obj = spatial.SpatialObject(array)

    # pylint guesses type wrong (self.obj.coords is a np.ndarray)
    # E1101: Instance of 'tuple' has no 'shape' member
    # pylint: disable=E1101
    def test_constructor(self):
        self.assertEqual(self.obj.coords.shape, (2, 3))
        assert_array_almost_equal(self.obj.coords[1], (1, 1, 1))

    def test_inertia_tensor(self):
        expected = np.full((3, 3), -0.5)
        np.fill_diagonal(expected, 1)
        assert_array_almost_equal(self.obj.inertia_tensor(), expected)


class TestRotation(unittest.TestCase):
    def test_rotation_X(self):
        # Rotation by 90°.
        ref = [
            [1.00, 0.00, 0.00, 0.00],
            [0.00, 0.00, -1.00, 0.00],
            [0.00, 1.00, 0.00, 0.00],
            [0.00, 0.00, 0.00, 1.00],
        ]
        m = linalg.rotation_matrix([90, 0, 0])
        assert_array_almost_equal(m, ref)

        # Rotation by 180°.
        ref = [
            [1.00, 0.00, 0.00, 0.00],
            [0.00, -1.00, -0.00, 0.00],
            [0.00, 0.00, -1.00, 0.00],
            [0.00, 0.00, 0.00, 1.00],
        ]
        m = linalg.rotation_matrix([180, 0, 0])
        assert_array_almost_equal(m, ref)

        # Rotation by 10°.
        ref = [
            [1.0000000, 0.000000, 0.000000, 0.000000],
            [0.0000000, 0.984808, -0.173648, 0.000000],
            [0.0000000, 0.173648, 0.984808, 0.000000],
            [0.0000000, 0.000000, 0.000000, 1.000000],
        ]
        m = linalg.rotation_matrix([10, 0, 0])
        assert_array_almost_equal(m, ref)

    def test_rotation_Y(self):
        # Rotation by 90°.
        ref = [
            [0.00, 0.00, 1.00, 0.00],
            [0.00, 1.00, 0.00, 0.00],
            [-1.00, 0.00, 0.00, 0.00],
            [0.00, 0.00, 0.00, 1.00],
        ]
        m = linalg.rotation_matrix([0, 90, 0])
        assert_array_almost_equal(m, ref)

        # Rotation by 180°.
        ref = [
            [-1.00, 0.00, 0.00, 0.00],
            [0.00, 1.00, 0.00, 0.00],
            [-0.00, 0.00, -1.00, 0.00],
            [0.00, 0.00, 0.00, 1.00],
        ]
        m = linalg.rotation_matrix([0, 180, 0])
        assert_array_almost_equal(m, ref)

        # Rotation by 10°.
        ref = [
            [0.9848078, 0.000000, 0.173648, 0.000000],
            [0.0000000, 1.000000, 0.000000, 0.000000],
            [-0.1736482, 0.000000, 0.984808, 0.000000],
            [0.0000000, 0.000000, 0.000000, 1.000000],
        ]
        m = linalg.rotation_matrix([0, 10, 0])
        assert_array_almost_equal(m, ref)

    def test_rotation_Z(self):
        # Rotation by 90°.
        ref = [
            [0.00, -1.00, 0.00, 0.00],
            [1.00, 0.00, 0.00, 0.00],
            [0.00, 0.00, 1.00, 0.00],
            [0.00, 0.00, 0.00, 1.00],
        ]
        m = linalg.rotation_matrix([0, 0, 90])
        assert_array_almost_equal(m, ref)

        # Rotation by 180°.
        ref = [
            [-1.00, -0.00, 0.00, 0.00],
            [0.00, -1.00, 0.00, 0.00],
            [0.00, 0.00, 1.00, 0.00],
            [0.00, 0.00, 0.00, 1.00],
        ]
        m = linalg.rotation_matrix([0, 0, 180])
        assert_array_almost_equal(m, ref)

        # Rotation by 10°.
        ref = [
            [0.9848078, -0.173648, 0.000000, 0.000000],
            [0.1736482, 0.984808, 0.000000, 0.000000],
            [0.0000000, 0.000000, 1.000000, 0.000000],
            [0.0000000, 0.000000, 0.000000, 1.000000],
        ]
        m = linalg.rotation_matrix([0, 0, 10])
        assert_array_almost_equal(m, ref)

    def test_rotation_XYZ(self):
        """Test rotation in X, Y and Z is equivalent to successive rotation
        in X, then Y, then Z."""
        m = linalg.rotation_matrix([10, 10, 10])
        ref = (
            linalg.rotation_matrix([10, 0, 0])
            .dot(linalg.rotation_matrix([0, 10, 0]))
            .dot(linalg.rotation_matrix([0, 0, 10]))
        )
        assert_array_almost_equal(m, ref)

    def test_rotate_by(self):
        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = np.array(list([i + 1.0, i + 11.0, i + 21.0] for i in range(10)))

        # Rotate by 12° along X-axis.
        transform.rotate_by(coords, [12, 0, 0])

        # Reference coordinates calculated with VMD::
        #     set m [transaxis x 12]
        #     set sel [atomselect 0 "all"]
        #     $sel move $m
        ref_coords = [
            [1.0000000, 6.393478, 22.828129],
            [2.0000000, 7.163715, 24.014187],
            [3.0000000, 7.933951, 25.200247],
            [4.0000000, 8.704186, 26.386307],
            [5.0000000, 9.474422, 27.572367],
            [6.0000000, 10.244658, 28.758427],
            [7.0000000, 11.014894, 29.944485],
            [8.0000000, 11.785130, 31.130545],
            [9.0000000, 12.555366, 32.316605],
            [10.0000000, 13.325602, 33.502663],
        ]

        assert_array_almost_equal(coords, ref_coords, decimal=5)

    def test_rotate_spatial_object(self):
        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = np.array(list([i + 1.0, i + 11.0, i + 21.0] for i in range(10)))

        # Initialize SpatialObject.
        o = spatial.SpatialObject(coords)

        # Rotate by 12° along X-axis.
        o.rotate_by([12, 0, 0])

        # Reference coordinates calculated with VMD::
        #     set m [transaxis x 12]
        #     set sel [atomselect 0 "all"]
        #     $sel move $m
        ref_coords = [
            [1.0000000, 6.393478, 22.828129],
            [2.0000000, 7.163715, 24.014187],
            [3.0000000, 7.933951, 25.200247],
            [4.0000000, 8.704186, 26.386307],
            [5.0000000, 9.474422, 27.572367],
            [6.0000000, 10.244658, 28.758427],
            [7.0000000, 11.014894, 29.944485],
            [8.0000000, 11.785130, 31.130545],
            [9.0000000, 12.555366, 32.316605],
            [10.0000000, 13.325602, 33.502663],
        ]

        assert_array_almost_equal(o.coords, ref_coords, decimal=5)

    def test_rotate_abstract_euler(self):
        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = np.array(list([i + 1.0, i + 11.0, i + 21.0] for i in range(10)))

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
        ref = np.array(
            [
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
            ]
        )

        assert_array_almost_equal(o.coords, ref)


class TestTransformation(unittest.TestCase):
    """Test transformations thanks to 4 x 4 matrices."""

    def test_transform_rotation(self):
        """Test transformation with a rotation matrix.

        In this test, the last row and columns are identity.
        """
        # Matrix calculated with VMD that correspond to a rotation
        # in the X-axis by 10°.
        m = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 0.984807753012208, -0.17364817766693033, 0.0],
                [0.0, 0.17364817766693033, 0.984807753012208, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )

        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = np.array(list([i + 1.0, i + 11.0, i + 21.0] for i in range(10)))

        ref = coords.copy()
        transform.rotate_by(ref, [10, 0, 0])

        transform.transform(coords, m)
        assert_array_almost_equal(coords, ref)

    def test_transform(self):
        """Test transformation with a random matrix."""
        # Input transformation matrix (just a random matrix).
        m = np.array(
            [
                [0.72898284, 0.43879558, 0.31011219, 0.32696626],
                [0.74291592, 0.00244003, 0.20881428, 0.91014385],
                [0.71045798, 0.85828462, 0.71570462, 0.96008097],
                [0.43315630, 0.30949433, 0.93486660, 0.48177513],
            ]
        )

        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = np.array(list([i + 1.0, i + 11.0, i + 21.0] for i in range(10)))

        # Reference calculated with VMD
        # with input topology = tests/data/test_10atoms.pdb.
        ref = np.array(
            [
                [0.5175050497055054, 0.253219336271286, 1.0914303064346313],
                [0.5412969589233398, 0.2738751769065857, 1.1091270446777344],
                [0.5621657967567444, 0.291993111371994, 1.1246496438980103],
                [0.5806189775466919, 0.3080138564109802, 1.1383754014968872],
                [0.5970528721809387, 0.3222815990447998, 1.1505991220474243],
                [0.611781895160675, 0.33506909012794495, 1.1615549325942993],
                [0.6250582337379456, 0.3465954065322876, 1.171429991722107],
                [0.6370866894721985, 0.35703831911087036, 1.1803770065307617],
                [0.6480352282524109, 0.3665436804294586, 1.1885206699371338],
                [0.6580431461334229, 0.37523239850997925, 1.1959648132324219],
            ]
        )

        transform.transform(coords, m)
        assert_array_almost_equal(coords, ref)

    def test_transform_spatial_object(self):
        """Test spatial.SpatialObject.transform."""
        # Input transformation matrix (just a random matrix).
        m = np.array(
            [
                [0.72898284, 0.43879558, 0.31011219, 0.32696626],
                [0.74291592, 0.00244003, 0.20881428, 0.91014385],
                [0.71045798, 0.85828462, 0.71570462, 0.96008097],
                [0.43315630, 0.30949433, 0.93486660, 0.48177513],
            ]
        )

        # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
        coords = np.array(list([i + 1.0, i + 11.0, i + 21.0] for i in range(10)))

        # Initialize SpatialObject.
        o = spatial.SpatialObject(coords)

        # Rotate by 12° along X-axis.
        o.transform(m)

        # Reference calculated with VMD
        # with input topology = tests/data/test_10atoms.pdb.
        ref = np.array(
            [
                [0.5175050497055054, 0.253219336271286, 1.0914303064346313],
                [0.5412969589233398, 0.2738751769065857, 1.1091270446777344],
                [0.5621657967567444, 0.291993111371994, 1.1246496438980103],
                [0.5806189775466919, 0.3080138564109802, 1.1383754014968872],
                [0.5970528721809387, 0.3222815990447998, 1.1505991220474243],
                [0.611781895160675, 0.33506909012794495, 1.1615549325942993],
                [0.6250582337379456, 0.3465954065322876, 1.171429991722107],
                [0.6370866894721985, 0.35703831911087036, 1.1803770065307617],
                [0.6480352282524109, 0.3665436804294586, 1.1885206699371338],
                [0.6580431461334229, 0.37523239850997925, 1.1959648132324219],
            ]
        )

        assert_array_almost_equal(o.coords, ref)


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
    def setUp(self):
        array = np.array(((0, 0, 0), (1, 1, 1)), dtype=float)
        self.obj = spatial.SpatialObject(array)

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


