# Python core libraries.
import math
import random

# Unit-test libraries.
import pytest
from pytest import approx

# Scientific libraries.
import numpy as np

# PTools imports.
from ptools import linalg
from ptools.linalg.matrix import (
    rotation_matrix,
    translation_matrix,
    transformation_matrix,
)

# More test-specific imports.
from .testing import assert_array_almost_equal


def generate_random_array(
    low: float = -10.0, high: float = 10.0, shape: int | tuple[int] = 1
) -> np.ndarray:
    """Generates a random array of floats."""
    return np.random.uniform(low, high, shape)


def generate_random_coordinates() -> np.ndarray:
    """Generate a coordinates-like array.

    The output array's size is random as well.
    """
    N = random.randint(2, 3000)
    return generate_random_array(-1000, 1000, (N, 3))


def generate_empty_coordinates() -> np.ndarray:
    """Returns an empty array."""
    return np.zeros((0, 3))


def test_angle():
    u = (2, 2, 0)
    v = (0, 3, 0)
    angle = linalg.angle(u, v)
    assert math.degrees(angle) == approx(45)


class TestDistanceToAxis:
    """Namespace that holds unit-tests for `ptools.linalg.distance_to_axis`."""

    def test_center_true(self):
        x = generate_random_coordinates()
        axis = generate_random_array(shape=(3,))
        actual = linalg.distance_to_axis(x, axis, center=True)
        expected = linalg.distance_to_axis(x - np.mean(axis, axis=0), axis, center=False)
        assert actual == approx(expected)

    def test_center_array(self):
        x = generate_random_coordinates()
        axis = generate_random_array(shape=(3,))
        center = generate_random_array(shape=(3,))
        actual = linalg.distance_to_axis(x, axis, center=center)
        expected = linalg.distance_to_axis(x - center, axis, center=False)
        assert actual == approx(expected)

    def test_center_false(self):
        x = generate_random_coordinates()
        axis = generate_random_array(shape=(3,))
        actual = linalg.distance_to_axis(x, axis, center=False)
        expected = np.linalg.norm(np.cross(x, axis))
        assert actual == approx(expected)


class TestTranslationMatrix:
    """Namespace that holds unit-tests for `ptools.linalg.translation_matrix`."""

    def test_vector(self):
        x = generate_random_array(shape=(3,))
        actual = translation_matrix(x)
        expected = np.array(
            [
                [1, 0, 0, x[0]],
                [0, 1, 0, x[1]],
                [0, 0, 1, x[2]],
                [0, 0, 0, 1],
            ]
        )
        assert_array_almost_equal(actual, expected)

    def test_scalar(self):
        x = random.random()
        actual = translation_matrix(x)
        expected = np.array(
            [
                [1, 0, 0, x],
                [0, 1, 0, x],
                [0, 0, 1, x],
                [0, 0, 0, 1],
            ]
        )
        assert_array_almost_equal(actual, expected)


class TestRotationMatrix:
    """Namespace that holds unit-tests for `ptools.linalg.rotation_matrix`."""

    def test_rotation_x_by_90(self):
        expected = [
            [1.00, 0.00, 0.00],
            [0.00, 0.00, -1.00],
            [0.00, 1.00, 0.00],
        ]
        actual = linalg.rotation_matrix([90, 0, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_x_by_180(self):
        expected = [
            [1.00, 0.00, 0.00],
            [0.00, -1.00, -0.00],
            [0.00, 0.00, -1.00],
        ]
        actual = linalg.rotation_matrix([180, 0, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_x_by_10(self):
        expected = [
            [1.0000000, 0.000000, 0.000000],
            [0.0000000, 0.984808, -0.173648],
            [0.0000000, 0.173648, 0.984808],
        ]
        actual = linalg.rotation_matrix([10, 0, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_y_by_90(self):
        expected = [
            [0.00, 0.00, 1.00],
            [0.00, 1.00, 0.00],
            [-1.00, 0.00, 0.00],
        ]
        actual = linalg.rotation_matrix([0, 90, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_y_by_180(self):
        expected = [
            [-1.00, 0.00, 0.00],
            [0.00, 1.00, 0.00],
            [-0.00, 0.00, -1.00],
        ]
        actual = linalg.rotation_matrix([0, 180, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_y_by_10(self):
        expected = [
            [0.9848078, 0.000000, 0.173648],
            [0.0000000, 1.000000, 0.000000],
            [-0.1736482, 0.000000, 0.984808],
        ]
        actual = linalg.rotation_matrix([0, 10, 0])
        assert_array_almost_equal(actual, expected)

    def test_rotation_z_by_90(self):
        expected = [
            [0.00, -1.00, 0.00],
            [1.00, 0.00, 0.00],
            [0.00, 0.00, 1.00],
        ]
        actual = linalg.rotation_matrix([0, 0, 90])
        assert_array_almost_equal(actual, expected)

    def test_rotation_z_by_180(self):
        expected = [
            [-1.00, -0.00, 0.00],
            [0.00, -1.00, 0.00],
            [0.00, 0.00, 1.00],
        ]
        actual = linalg.rotation_matrix([0, 0, 180])
        assert_array_almost_equal(actual, expected)

    def test_rotation_z_by_10(self):
        expected = [
            [0.9848078, -0.173648, 0.000000],
            [0.1736482, 0.984808, 0.000000],
            [0.0000000, 0.000000, 1.000000],
        ]
        actual = linalg.rotation_matrix([0, 0, 10])
        assert_array_almost_equal(actual, expected)

    def test_rotation_xyz(self):
        """Tests rotation order."""
        angles = [10, 50, 80]
        by_x = linalg.rotation_matrix([angles[0], 0, 0])
        by_y = linalg.rotation_matrix([0, angles[1], 0])
        by_z = linalg.rotation_matrix([0, 0, angles[2]])

        actual = linalg.rotation_matrix(angles)
        expected = by_z.dot(by_y).dot(by_x)
        assert_array_almost_equal(actual, expected)

    def test_rotation_zyx(self):
        """Tests rotation order."""
        angles = [10, 10, 10]
        by_x = linalg.rotation_matrix([angles[0], 0, 0])
        by_y = linalg.rotation_matrix([0, angles[1], 0])
        by_z = linalg.rotation_matrix([0, 0, angles[2]])

        actual = linalg.rotation_matrix(angles, sequence="zyx")
        expected = by_x.dot(by_y).dot(by_z)
        assert_array_almost_equal(actual, expected)


class TestTransformationMatrixFromVectors:
    def test_rotations_only(self):
        rotation = generate_random_array(shape=(3,))

        expected = np.eye(4)
        expected[:3, :3] = rotation_matrix(rotation)

        actual = transformation_matrix(rotation=rotation)
        assert actual.shape == (4, 4)
        assert_array_almost_equal(actual, expected)

    def test_translations_only(self):
        translation = generate_random_array(shape=(3,))

        expected = np.eye(4)
        expected[:3, 3] = translation

        actual = transformation_matrix(translation=translation)
        assert actual.shape == (4, 4)
        assert_array_almost_equal(actual, expected)

    def test_rotation_and_translations(self):
        rotation = generate_random_array(shape=(3,))
        translation = generate_random_array(shape=(3,))

        expected = np.eye(4)
        expected[:3, :3] = rotation_matrix(rotation)
        expected[:3, 3] = translation

        actual = transformation_matrix(translation, rotation)
        assert actual.shape == (4, 4)
        assert_array_almost_equal(actual, expected)


class TestTransformationMatrixFromMatrices:
    def test_rotations_only(self):
        rotation = generate_random_array(shape=(3,))

        expected = np.eye(4)
        expected[:3, :3] = rotation_matrix(rotation)

        actual = transformation_matrix(rotation=rotation_matrix(rotation))
        assert actual.shape == (4, 4)
        assert_array_almost_equal(actual, expected)

    def test_translations_only(self):
        translation = generate_random_array(shape=(3,))

        expected = np.eye(4)
        expected[:3, 3] = translation

        actual = transformation_matrix(translation=translation_matrix(translation))
        assert actual.shape == (4, 4)
        assert_array_almost_equal(actual, expected)

    def test_rotation_and_translations(self):
        rotation = generate_random_array(shape=(3,))
        translation = generate_random_array(shape=(3,))

        expected = np.eye(4)
        expected[:3, :3] = rotation_matrix(rotation)
        expected[:3, 3] = translation

        actual = transformation_matrix(
            translation_matrix(translation), rotation_matrix(rotation)
        )
        assert actual.shape == (4, 4)
        assert_array_almost_equal(actual, expected)


# ======================================================================================
class TestRotationMatrixAroundAxis:
    def test_rotation_matrix_round_axis(self):
        # rotation of 90° along Z-axis.
        actual = linalg.rotation_matrix_around_axis([0, 0, 1], 90)
        expected = np.eye(4)
        expected[:3, :3] = rotation_matrix([0, 0, 90])
        assert_array_almost_equal(actual, expected)
        print(actual)
        # assert 1 == 2


# ======================================================================================
class TestAttractEulerRotationMatrix:
    def test_rotation_x_by_90(self):
        alpha = 90
        actual = linalg.matrix.attract_euler_rotation_matrix([alpha, 0, 0])

        # Attract convention: first angle is actually a rotation along the Z-axis.
        expected = rotation_matrix([0, 0, alpha], degrees=False)
        assert_array_almost_equal(actual, expected)

    def test_rotation_y_by_90(self):
        alpha = 90
        actual = linalg.matrix.attract_euler_rotation_matrix([0, alpha, 0])
        expected = rotation_matrix([0, alpha, 0], degrees=False)
        assert_array_almost_equal(actual, expected)

    def test_rotation_z_by_90(self):
        alpha = 90
        actual = linalg.matrix.attract_euler_rotation_matrix([0, 0, alpha])
        # Attract convention: matrix has to be transposed... don't know what's going on here...
        expected = rotation_matrix([0, 0, alpha], degrees=False).T
        assert_array_almost_equal(actual, expected)

    def test_rotation_xyz(self):
        # actually a non regression test
        angles = [10, 50, 80]
        actual = linalg.matrix.attract_euler_rotation_matrix(angles)
        expected = np.array(
            [
                [0.630074, 0.744674, 0.220151],
                [-0.775995, 0.614376, 0.142737],
                [-0.028963, -0.260771, 0.964966],
            ]
        )
        assert_array_almost_equal(actual, expected)

    def test_legacy(self):
        angles = [10, 50, 80]
        actual = linalg.matrix.attract_euler_rotation_matrix_legacy(*angles)
        expected = linalg.matrix.attract_euler_rotation_matrix(angles)
        assert_array_almost_equal(actual, expected)


# ======================================================================================
class TestCentroid:
    """Namespace that holds unit-tests for ptools.linalg.centroid."""

    def test_random_coordinates(self):
        x = generate_random_coordinates()
        assert_array_almost_equal(linalg.centroid(x), np.mean(x, axis=0))

    def test_empty_coordinates(self):
        x = generate_empty_coordinates()
        err = "cannot compute the centroid of an empty array"
        with pytest.raises(ZeroDivisionError, match=err):
            linalg.centroid(x)


# ======================================================================================
class TestCenterOfMass:
    """Namespace that holds unit-tests for ptools.linalg.center_of_mass."""

    def test_random_coordinates(self):
        x = generate_random_coordinates()
        w = generate_random_array(shape=x.shape[0])
        actual = linalg.center_of_mass(x, w)
        expected = np.average(x, axis=0, weights=w)
        assert_array_almost_equal(actual, expected)

    def test_empty_coordinates(self):
        x = generate_empty_coordinates()
        w = np.ones(x.shape[0])
        err = "cannot compute the center of mass of an empty array"
        with pytest.raises(ZeroDivisionError, match=err):
            linalg.center_of_mass(x, w)

    def test_empty_weights(self):
        x = generate_random_coordinates()
        w = np.ones(0)
        err = "input array and weights should be the same size"
        with pytest.raises(ValueError, match=err):
            linalg.center_of_mass(x, w)


# ======================================================================================

# pylint: disable=R0903
class TestInertiaTensor:
    """Namespace that holds unit-tests for ptools.linalg.inertia_tensor."""

    def test_fails_if_array_is_not_n_by_3(self):
        x = np.array(((0.0, 0.0), (1.0, 1.0)))
        w = np.ones((x.shape[0]))
        err = "inertia tensor can only be calculated on a N x 3 array"
        with pytest.raises(ValueError, match=err):
            linalg.inertia_tensor(x, w)


class TestDNA:
    def test_adna(self):
        # non-regression test
        expected = np.array(
            [
                [0.85554587, -0.51531864, 0.04987943, 0.9404263],
                [0.51647696, 0.85620086, -0.01310093, -2.39183975],
                [-0.03595566, 0.03697002, 0.99866932, 3.30599325],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        actual = linalg.matrix.adna_tranformation_matrix()
        assert_array_almost_equal(actual, expected)

    def test_bdna(self):
        # non-regression test
        expected = np.array(
            [
                [8.09100476e-01, -5.84986445e-01, -5.61006103e-02, 2.59187817e-01],
                [5.85826775e-01, 8.10434321e-01, -1.78910064e-03, -1.34486756e00],
                [4.65124597e-02, -3.14176774e-02, 9.98423518e-01, 3.28877603e00],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        actual = linalg.matrix.bdna_tranformation_matrix()
        assert_array_almost_equal(actual, expected)
