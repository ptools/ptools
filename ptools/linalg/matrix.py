# Python core libraries.
import math

# Scientific libraries.
import numpy as np
from numpy.typing import ArrayLike
import scipy.linalg as L


def translation_matrix(direction: ArrayLike = np.zeros(3)) -> np.ndarray:
    """ReturnS the matrix to translate by direction vector."""
    matrix = np.identity(4)
    if np.isscalar(direction):
        matrix[:3, 3] = direction
    else:
        matrix[:3, 3] = direction[:3]
    return matrix


def transformation_matrix(
    translation: ArrayLike = np.zeros(3), rotation: ArrayLike = np.zeros(3)
) -> np.ndarray:
    """Returns a 4x4 transformation matrix."""
    matrix = rotation_matrix(rotation)
    matrix[:3, 3] = translation
    return matrix


def rotation_matrix(angles: np.ndarray = np.zeros(3)) -> np.ndarray:
    """Return the rotation matrix around the X, Y and Z axes.

    The matrix rotates first along X-axis, then Y, then Z.

    Args:
        angles ([float, float, float]): rotation angle around
                the X-, Y- and Z-axes

    Returns:
        numpy.ndarray: 4x4 rotation matrix.
    """
    assert np.shape(angles) == (3,)

    alpha = math.radians(angles[0])
    beta = math.radians(angles[1])
    gamma = math.radians(angles[2])
    r = np.identity(4)
    r[0, 0] = math.cos(beta) * math.cos(gamma)
    r[0, 1] = -math.cos(beta) * math.sin(gamma)
    r[0, 2] = math.sin(beta)

    r[1, 0] = math.cos(alpha) * math.sin(gamma) + math.sin(alpha) * math.sin(
        beta
    ) * math.cos(gamma)
    r[1, 1] = math.cos(alpha) * math.cos(gamma) - math.sin(alpha) * math.sin(
        beta
    ) * math.sin(gamma)
    r[1, 2] = -math.sin(alpha) * math.cos(beta)

    r[2, 0] = math.sin(alpha) * math.sin(gamma) - math.cos(alpha) * math.sin(
        beta
    ) * math.cos(gamma)
    r[2, 1] = math.sin(alpha) * math.cos(gamma) + math.cos(alpha) * math.sin(
        beta
    ) * math.sin(gamma)
    r[2, 2] = math.cos(alpha) * math.cos(beta)
    return r


def rotation_matrix_around_axis(
    axis: np.ndarray, amount: float, center: np.ndarray = np.zeros(3)
):
    """Returns the rotation matrix to rotate around the axis by given angle.

    Args:
        axis (np.ndarray (N x 3)): axis of rotation
        amount (float): amount in radians
        center (np.ndarray (3 x 1)): center of rotation

    Returns:
        np.ndarray: 4 x 4 rotation matrix
    """
    assert np.shape(axis) == (3,)
    assert np.shape(center) == (3,)

    origin_matrix = translation_matrix(-center)
    offset_matrix = translation_matrix(+center)
    rotation = L.expm(np.cross(np.identity(3), axis / L.norm(axis) * amount))
    matrix = np.identity(4)
    matrix[:3, :3] = rotation
    return np.matmul(np.matmul(offset_matrix, matrix), origin_matrix)


def attract_euler_rotation_matrix(phi: float, ssi: float, rot: float):
    """Return the rotation matrix for an Euler rotation with Attract
    convention.

    Args:
        phi (float):
        ssi (float):
        rot (float):

    Returns:
        numpy.ndarray: 3 x 3 matrix
    """
    cs = math.cos(ssi)
    cp = math.cos(phi)
    ss = math.sin(ssi)
    sp = math.sin(phi)
    crot = math.cos(rot)
    srot = math.sin(rot)

    cscp = cs * cp
    cssp = cs * sp
    sscp = ss * cp
    sssp = ss * sp

    matrix = np.identity(3)

    matrix[0][0] = crot * cscp + srot * sp
    matrix[0][1] = srot * cscp - crot * sp
    matrix[0][2] = sscp

    matrix[1][0] = crot * cssp - srot * cp
    matrix[1][1] = srot * cssp + crot * cp
    matrix[1][2] = sssp

    matrix[2][0] = -crot * ss
    matrix[2][1] = -srot * ss
    matrix[2][2] = cs

    return matrix


def orientation_matrix(
    coords: np.ndarray, vector: np.ndarray, target: np.ndarray
) -> np.ndarray:
    """Calculates an orientation matrix.

    Orients vector on target.

    IMPORTANT: does not take into account the center of mass.

    Args:
        coords (np.ndarray <N x 3>): coordinates
        vector (np.ndarray (1 x 3)): reference for rotation
        target (np.ndarray (1 x 3)): target for rotation

    Returns:
        np.ndarray (4 x 4): transformation matrix

    Example:
        >>> # Align the 3rd principal axis of inertia on the Z-axis.
        >>> I = inertia_tensors(rigidbody)
        >>> T = orientation_matrix(rigidbody.coords, I[2], [0, 0, 1])
        >>> receptor.transform(T)
    """
    assert np.shape(vector) == np.shape(target)
    assert np.ndim(coords) == 2
    assert np.ndim(vector) == 1 and len(vector) == np.shape(coords)[1]

    com = coords.mean(axis=0)

    vec1 = vector / np.linalg.norm(vector)
    vec2 = target / np.linalg.norm(target)

    axis = np.cross(vec1, vec2)
    sine = np.linalg.norm(axis)
    cosine = np.dot(vec1, vec2)
    amount = math.atan2(sine, cosine)

    return rotation_matrix_around_axis(axis, amount, center=com)


def ab_rotation_matrix(A: np.ndarray, B: np.ndarray, amount: float) -> np.ndarray:
    """Returns the rotation matrix to rotate around axis (A, B) by amount (in radians)."""
    return rotation_matrix_around_axis(B - A, amount, A)
