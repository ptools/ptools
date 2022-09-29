"""Move 3D-objects."""

# Scientific libraries.
import numpy as np
from numpy.typing import ArrayLike

# PTools libraries
from .matrix import (
    attract_euler_rotation_matrix,
    orientation_matrix,
    rotation_matrix,
    rotation_matrix_around_axis,
)


def translate(coords: np.ndarray, t: ArrayLike):
    """In-place translation of coordinates."""

    def _translate_scalar(x):
        np.add(coords, x, coords)

    def isscalar(s):
        return isinstance(s, (float, int))

    if isscalar(t):
        _translate_scalar(t)
    else:
        t = np.array(t)
        if t.shape == (4, 4):
            t = t[:3, 3]
        elif t.shape != (3,):
            raise ValueError(
                "Dimensions error: expected 3 x 1 or 4 x 4 " f"(got {t.shape})"
            )
        np.add(coords, t, coords)


def rotate_by(coords: np.ndarray, angles: ArrayLike = np.zeros(3)):
    """In-place rotation of coordinates around X, Y and Z axes.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        angles ([float, float, float]): rotation angle around
                the X-, Y- and Z-axes
    """
    matrix = rotation_matrix(angles)
    rotate(coords, matrix)


def rotate(coords: np.ndarray, r: ArrayLike):
    """In-place rotation of coordinates using a rotation matrix or 3 angles."""
    r = np.array(r)
    if r.shape == (3,):
        matrix = rotation_matrix(r)
    elif r.shape in ((4, 4), (3, 3)):
        matrix = r
    else:
        raise ValueError(
            "Dimensions error: expected 3 x 1 or 4 x 4 " f"(got {r.shape}) "
        )

    coords[:] = np.inner(coords, matrix[:3, :3])


def transform(coords: np.ndarray, matrix: ArrayLike):
    """In-place transformation of coordinates by a 4 x 4 matrix.

    This function should be used only when the transformation matrix
    last row and last column differ from (0, 0, 0, 1).
    When its not the case, just convert the transformation matrix to 3 x 3
    and multiply coordinates::

        >>> coords[:] = numpy.inner(coords, matrix[:3, :3])

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        matrix (numpy.ndarray): 4 x 4 matrix
    """
    # Reshape coords to N x 4 with the last component being 1.
    n = coords.shape[1]
    a = np.ones((coords.shape[0], n + 1))
    a[:, :-1] = coords

    # Apply transformation matrix to coordinates.
    a[:] = np.inner(a, matrix)

    # Weight coordinates and reshape into N x 3 again.
    coords[:] = (a[:, 0:n].T / a[:, n]).T


def orient(coords: np.ndarray, vector: ArrayLike, target: ArrayLike):
    """Orients coordinates.

    Args:
        coords (np.ndarray (N x 3)): coordinate array:
        vector (np.ndarray (3 x 1)): source vector
        target (np.ndarray (3 x 1)): target vector
    """
    matrix = orientation_matrix(coords, vector, target)
    transform(coords, matrix)


def attract_euler_rotate(coords: np.ndarray, angles: ArrayLike = np.zeros(3)):
    """In-place Euler rotation of coordinates using the Attract convention."""
    matrix = attract_euler_rotation_matrix(angles)
    coords[:] = np.inner(coords, matrix)


def ab_rotate(coords: np.ndarray, A: ArrayLike, B: ArrayLike, amount: float, degrees: bool = True):
    """Rotates coords around axis (A, B) by amount theta (in radians)."""
    matrix = rotation_matrix_around_axis(B - A, amount, A, degrees)
    transform(coords, matrix)
