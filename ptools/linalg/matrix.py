# Python core libraries.
import math

# Scientific libraries.
import numpy as np
import scipy

# Type hinting libraries.
from .._typing import ArrayLike


def translation_matrix(direction: ArrayLike = np.zeros(3)) -> np.ndarray:
    """Returns a 4 x 4 matrix to translate by direction vector."""
    matrix = np.identity(4)
    if np.isscalar(direction):
        matrix[:3, 3] = direction
    else:
        matrix[:3, 3] = direction[:3]
    return matrix


def transformation_matrix(
    translation: ArrayLike = np.zeros(3), rotation: ArrayLike = np.zeros(3)
) -> np.ndarray:
    """Returns a 4x4 transformation matrix.

    Args:
        translations: offsets for X- Y- Z-axes translations or 4 x 4 translation matrix
        rotation: angles for X- Y- Z-axes rotations or 3 x 3 rotation matrix
    """
    matrix = np.identity(4)
    translation = np.asarray(translation)
    rotation = np.asarray(rotation)

    if np.shape(rotation) == (3,):
        matrix[:3, :3] = rotation_matrix(rotation)
    elif np.shape(rotation) == (3, 3):
        matrix[:3, :3] = rotation
    else:
        raise ValueError(
            f"rotation expects either vector of size 3 or 3 x 3 rotation matrix (found {np.shape(rotation)})"
        )

    if np.shape(translation) == (3,):
        matrix[:3, 3] = translation
    elif np.shape(translation) == (4, 4):
        matrix[:, 3] = translation[:, 3]
    else:
        raise ValueError(
            f"translation expects either vector of size 3 or 3 x 3 translation matrix (found {np.shape(translation)})"
        )
    return matrix


def rotation_matrix(
    angles: ArrayLike = np.zeros(3), sequence: str = "xyz", degrees: bool = True
) -> np.ndarray:
    """Returns a rotation matrix in from Euler angles.

    This is an alias for `scipy.spatial.transform.Rotation.from_euler`.

    Args:
        angles ([float, float, float]): rotation angle around
                the X-, Y- and Z-axes
        sequence: order in which rotations are applied
        degrees: angles are given either in degrees or in radians

    Returns:
        numpy.ndarray: 3 x 3 rotation matrix.
    """
    return scipy.spatial.transform.Rotation.from_euler(
        sequence, angles, degrees
    ).as_matrix()


def rotation_matrix_around_axis(
    axis: np.ndarray, amount: float, center: np.ndarray = np.zeros(3), degrees=True
):
    """Returns the rotation matrix to rotate around the axis by given angle.

    Args:
        axis (3 x 1): axis of rotation
        amount: angle
        center (3 x 1): center of rotation
        degrees: amount is given in degrees (if False, it is given in radians)

    Returns:
        np.ndarray: 4 x 4 transformation matrix
    """
    assert np.shape(axis) == (3,)
    assert np.shape(center) == (3,)

    if degrees:
        amount = math.radians(amount)

    origin_matrix = translation_matrix(-np.array(center))
    offset_matrix = translation_matrix(+np.array(center))
    rotation = scipy.linalg.expm(
        np.cross(np.identity(3), axis / scipy.linalg.norm(axis) * amount)
    )
    matrix = np.identity(4)
    matrix[:3, :3] = rotation
    return offset_matrix @ matrix @ origin_matrix


def attract_euler_rotation_matrix(angles: ArrayLike = np.zeros(3)) -> np.ndarray:
    """Return the rotation matrix for an Euler rotation with Attract
    convention.

    Args:
        angles: 3 angles for rotation on Z-, Y- and X- axes, respectively, given in radians.

    Returns:
        numpy.ndarray: 3 x 3 matrix
    """
    by_x = rotation_matrix(
        [0, 0, angles[2]], degrees=False
    ).T  # what the hell is going on here???
    by_y = rotation_matrix([0, angles[1], 0], degrees=False)
    by_z = rotation_matrix([0, 0, angles[0]], degrees=False)
    return by_z @ by_y @ by_x


def attract_euler_rotation_matrix_legacy(phi: float, ssi: float, rot: float):
    """Legagy version of `attract_euler_rotation_matrix`.

    This function will be removed in the future.

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
    coords: np.ndarray, vector: ArrayLike, target: ArrayLike
) -> np.ndarray:
    """Calculates an orientation matrix.

    Orients vector on target.

    IMPORTANT: does not take into account the center of mass.

    Args:
        coords (np.ndarray <N x 3>): coordinates
        vector (ArrayLike (1 x 3)): reference for rotation
        target (ArrayLike (1 x 3)): target for rotation

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

    vec1 = np.asarray(vector) / np.linalg.norm(vector)
    vec2 = np.asarray(target) / np.linalg.norm(target)

    axis = np.cross(vec1, vec2)
    sine = np.linalg.norm(axis)
    cosine = np.dot(vec1, vec2)
    amount = math.atan2(sine, cosine)

    return rotation_matrix_around_axis(axis, amount, center=com, degrees=False)


def ab_rotation_matrix(A: np.ndarray, B: np.ndarray, amount: float, degrees=True) -> np.ndarray:
    """Returns the rotation matrix to rotate around axis (A, B) by amount."""
    return rotation_matrix_around_axis(B - A, amount, A, degrees)


def shift_matrix(offset: float) -> np.ndarray:
    """Returns a shift matrix, which is a translation matrix along the X-axis."""
    return translation_matrix([offset, 0, 0])


def slide_matrix(offset: float) -> np.ndarray:
    """Returns a slide matrix, which is a translation matrix along the Y-axis."""
    return translation_matrix([0, offset, 0])


def rise_matrix(offset: float) -> np.ndarray:
    """Returns a rise matrix, which is a translation matrix along the Z-axis."""
    return translation_matrix([0, 0, offset])


def tilt_matrix(alpha: float) -> np.ndarray:
    """Returns a tilt matrix, which is a rotation matrix along the X-axis."""
    return transformation_matrix(rotation=rotation_matrix([alpha, 0, 0]))


def roll_matrix(alpha: float) -> np.ndarray:
    """Returns a roll matrix, which is a rotation matrix along the Y-axis."""
    return transformation_matrix(rotation=rotation_matrix([0, alpha, 0]))


def twist_matrix(alpha: float) -> np.ndarray:
    """Returns a twist matrix, which is a rotation matrix along the Z-axis."""
    return transformation_matrix(rotation=rotation_matrix([0, 0, alpha]))


# =======================================================================================
# Matrices to build "classic" DNA.
#
# Parameters taken from
# S.Arnott, R.Chandrasekaran, D.L.Birdsall, A.G.W.Leslie and R.L.Ratliff, Nature 283, 743-746 (1980)


def adna_tranformation_matrix() -> np.ndarray:
    matrix = (
        twist_matrix(31.1185909091)
        @ roll_matrix(2.06055181818)
        @ tilt_matrix(2.12008054545)
        @ rise_matrix(3.37983727273)
        @ slide_matrix(-2.41029181818)
        @ shift_matrix(-0.549621454545)
    )
    return matrix


def bdna_tranformation_matrix() -> np.ndarray:
    matrix = (
        twist_matrix(35.9063052632)
        @ roll_matrix(-2.66592947368)
        @ tilt_matrix(-1.80234789474)
        @ rise_matrix(3.27145684211)
        @ slide_matrix(-1.34487389474)
        @ shift_matrix(-0.425181378947)
    )
    return matrix
