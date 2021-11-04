"""ptools.spatial - Defines classes and functions to work on spatial
(mostly 3D) data."""

# Allows to use type hinting of class within itself
# e.g. SpatialObject.dist(self, other: SpatialObject)
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import expm, norm


class SpatialObject:
    """An object with coordinates.

    Implements basic spatial operations such as translation, rotation, etc.
    """

    def __init__(self, coords: np.ndarray = np.zeros(3)):
        self._coords = coord3d(coords)

    def dist(self, other: SpatialObject) -> float:
        """Returns the euclidean distance between two objects."""
        return dist(self.coords, other.coords)

    def copy(self) -> SpatialObject:
        """Returns a copy of itself."""
        return self.__class__(self.coords)

    @property
    def coords(self) -> np.ndarray:
        """Get SpatialObject cartesian coordinates."""
        return self._coords

    @coords.setter
    def coords(self, pos: np.ndarray | list[float]):
        """Set SpatialObject cartesian coordinates."""
        self._coords = coord3d(pos)

    def centroid(self) -> np.ndarray:
        """Returns an spatial object geometric center."""
        return centroid(self._coords)

    def translate(self, direction: int | float | np.ndarray):
        """Translate object coordinates using vector `direction`.

        Args:
            direction ((int, float) or np.ndarray): scalar, 1 x 3 shaped vector or
                4 x 4 matrix
        """
        translate(self.coords, direction)

    def center(self, origin: np.ndarray = np.zeros(3)):
        """Center spatial object on `origin`."""
        self.translate(origin - self.coords)

    def rotate_by(self, angles: np.ndarray):
        """Rotate object coordinates around X, Y and Z axes.

        Args:
            angles ([float, float, float]): rotation angle around
                the X-, Y- and Z-axes
        """
        rotate_by(self.coords, angles)

    def rotate(self, rotation: np.ndarray):
        """Rotate object using rotation matrix or 3 angles.

        Args:
            rotation (np.ndarray): 3 x 3 matrix of 3 x 1 vector of angles.
        """
        rotate(self.coords, rotation)

    def ab_rotate(self, A: np.ndarray, B: np.ndarray, amount: float):
        """PTools rotation around axis."""
        ab_rotate(self.coords, A, B, amount)

    def transform(self, matrix: np.ndarray):
        """Transform object by 4x4 matrix."""
        transform(self.coords, matrix)

    def move(self, matrix: np.ndarray):
        """Transform object by 4x4 matrix.

        This is an alias for `SpatialObject.transform`.
        """
        self.transform(matrix)

    def moveby(self, direction: np.ndarray | int):
        """Translate object coordinates using vector or scalar `direction`.

        This is an alias for `SpatialObject.translate`.
        """
        self.translate(direction)

    def attract_euler_rotate(self, phi: float, ssi: float, rot: float):
        """Rotate object with Attract convention."""
        attract_euler_rotate(self.coords, phi, ssi, rot)

    def orient(
        self, vector: np.ndarray | list[float], target: np.ndarray | list[float]
    ):
        """Orients a SpatialObject."""
        orient(self.coords, vector, target)

    def tensor_of_inertia(
        self, weights: np.ndarray = None, method: str = "fast"
    ) -> np.ndarray:
        """Returns a SpatialObject inertial tensors."""
        return tensor_of_inertia(self.coords, weights, method)

    def distance_to_axis(self, axis: np.ndarray) -> float:
        """Returns the SpatialObject distance to an arbitrary axis."""
        return np.linalg.norm(np.cross(self.coords, axis))


# pylint: disable-msg=W1113
def coord3d(value: np.ndarray = np.zeros(3), *args) -> np.ndarray:
    """Convert an iterable of size 3 to a 1 x 3 shaped numpy array of floats.

    Can be called either with an iterable, either with 3 arguments.

    Examples:
        >>> coord3d((1, 2, 3))
        array([1., 2., 3.])
        >>> coord3d(1, 2, 3)
        array([1., 2., 3.])
        >>> coord3d(2)
        array([2., 2., 2.])
    """
    # Checks function was called with adequate number of arguments.
    if args:
        if len(args) != 2:
            raise TypeError(
                f"Coordinates must be initialized either "
                f"with 1 or 3 arguments (found {value=}, {args=})"
            )

    # Initialization from a single value.
    if not args:
        # value is a scalar.
        if isinstance(value, (int, float)):
            return _coordinates_from_scalar(value)
        if isinstance(value, (list, tuple, np.ndarray)):
            return _coordinates_from_vector(value)
        # value is some invalid type.
        raise TypeError(f"Invalid coordinates initialization from {value}")

    # Initialization from 3 values.
    return np.array((value, *args))


def _coordinates_from_scalar(value: int | float) -> np.ndarray:
    """Returns a numpy array of size 3 filled with value."""
    return np.full((3,), value, dtype=float)


def _coordinates_from_vector(value: list[float] | np.ndarray) -> np.ndarray:
    """Returns a numpy array of size 3 initialized with a vector."""
    # value is an iterable: numpy will raise a ValueError if some element
    # is not compatible with float.
    if isinstance(value, (list, tuple)):
        array = np.array(value, dtype=float)
    elif isinstance(value, np.ndarray):
        array = value.astype("float64")
    # Checks array corresponds to 3D coordinates
    if len(array.shape) > 2:
        raise ValueError(
            "3D coordinate array should have at "
            f"most 2 dimensions (found {array.shape}"
        )
    if (len(array.shape) == 2 and array.shape[1] != 3) or (
        len(array.shape) == 1 and array.shape[0] != 3
    ):
        raise ValueError(
            "3D coordinate array should be N x 3 " f"(found {array.shape})"
        )
    return array


def centroid(x: np.ndarray) -> np.ndarray:
    """Return x centroid.

    Centroid is the average coordinate along axis 0.
    """
    return np.mean(x, axis=0)


def center_of_mass(coords: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Return the center of mass (barycenter)."""
    weights = weights.reshape(-1, 1)
    return (coords * weights).sum(axis=0) / weights.sum()


def angle(u: np.ndarray, v: np.ndarray) -> float:
    """Returns the angle between two vectors in radians."""
    return math.acos(np.dot(u, v) / (norm(u) * norm(v)))


def _tensor_of_inertia_accurate(coords: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Returns the inertia tensors of a set of atoms."""
    com = center_of_mass(coords, weights)
    X = coords - com

    Ixx = np.sum(weights * (X[:, 1] ** 2 + X[:, 2] ** 2))
    Iyy = np.sum(weights * (X[:, 0] ** 2 + X[:, 2] ** 2))
    Izz = np.sum(weights * (X[:, 0] ** 2 + X[:, 1] ** 2))
    Ixy = -np.sum(weights * X[:, 0] * X[:, 1])
    Iyz = -np.sum(weights * X[:, 1] * X[:, 2])
    Ixz = -np.sum(weights * X[:, 0] * X[:, 2])

    I = np.array(
        [
            [Ixx, Ixy, Ixz],
            [Ixy, Iyy, Iyz],
            [Ixz, Iyz, Izz],
        ]
    )

    return I


def _tensor_of_inertia_fast(coords: np.ndarray) -> np.ndarray:
    """Returns the inertia tensors of a set of atoms.

    Fast: does not take into account the weights
    """
    com = centroid(coords)
    X = coords - com
    tensors = np.dot(X.transpose(), X)
    return tensors


def tensor_of_inertia(
    coords: np.ndarray, weights: np.ndarray = None, method: str = "accurate"
) -> np.ndarray:
    """Returns the inertia tensors of a set of atoms.

    Allows to choose between "accurate" (weighted) and "fast" methods.
    """
    if method not in ("accurate", "fast"):
        raise ValueError(
            "parameter 'method' should be 'accurate' or 'fast' " f"(found {method=!r})"
        )
    if method == "accurate":
        if weights is None:
            raise ValueError("need weights to compute accurate tensor of " "inertia")
        return _tensor_of_inertia_accurate(coords, weights)
    return _tensor_of_inertia_fast(coords)  # method is "fast"


def principal_axes(tensor: np.ndarray, sort: bool = True) -> np.ndarray:
    """Return the principal axes given the tensor of inertia.

    Computes the eigenvalues and right eigenvectors of the tensor of inertia.

    Args:
        tensor (np.ndarray): 3 x 3 matrix
        sort (bool): sort axes by amplitude

    Returns:
        np.ndarray: 3 x 3 matrix
    """
    evalues, evectors = np.linalg.eig(tensor)
    if sort:
        order = np.argsort(evalues)[::-1]  # descending order
        evectors = evectors[:, order].transpose()
    return evectors


#
# Translation routines
#
def translation_matrix(direction: np.ndarray = np.zeros(3)) -> np.ndarray:
    """Return the matrix to translate by direction vector."""
    matrix = np.identity(4)
    matrix[:3, 3] = direction[:3]
    return matrix


def translate(coords: np.ndarray, t: np.ndarray | float | int):
    """In-place translation of coordinates by a vector."""

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


#
#  Rotation routines
#
def rotation_matrix(angles: np.ndarray = np.zeros(3)) -> np.ndarray:
    """Return the rotation matrix around the X, Y and Z axes.

    The matrix rotates first along X-axis, then Y, then Z.

    Args:
        angles ([float, float, float]): rotation angle around
                the X-, Y- and Z-axes

    Returns:
        numpy.ndarray: 4x4 rotation matrix.
    """
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


def rotate_by(coords: np.ndarray, angles: np.ndarray = np.zeros(3)):
    """In-place rotation of coordinates around X, Y and Z axes.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        angles ([float, float, float]): rotation angle around
                the X-, Y- and Z-axes
    """
    matrix = rotation_matrix(angles)
    rotate(coords, matrix)


def rotate(coords: np.ndarray, r: np.ndarray):
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


def attract_euler_rotation_matrix(phi: float, ssi: float, rot: float):
    """Return the rotation matrix for an Euler rotation with Attract
    convention.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
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


def attract_euler_rotate(coords: np.ndarray, phi: float, ssi: float, rot: float):
    """In-place Euler rotation of coordinates with the Attract convention.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        phi (float):
        ssi (float):
        rot (float):
    """
    matrix = attract_euler_rotation_matrix(phi, ssi, rot)
    coords[:] = np.inner(coords, matrix)


def ab_rotation_matrix(A: np.ndarray, B: np.ndarray, amount: float) -> np.ndarray:
    """Returns the rotation matrix to rotate around axis (A, B) by amount (in radians)."""
    return rotation_matrix_around_axis(B - A, amount, A)


def ab_rotate(coords: np.ndarray, A: np.ndarray, B: np.ndarray, amount: float):
    """Rotates coords around axis (A, B) by amount theta (in radians)."""
    matrix = rotation_matrix_around_axis(B - A, amount, A)
    transform(coords, matrix)


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

    def _translation_matrix(direction):
        matrix = np.identity(4)
        matrix[:3, 3] = direction[:3]
        return matrix

    origin_matrix = _translation_matrix(-center)
    offset_matrix = _translation_matrix(+center)
    rotation = expm(np.cross(np.identity(3), axis / norm(axis) * amount))
    matrix = np.identity(4)
    matrix[:3, :3] = rotation
    return np.matmul(np.matmul(offset_matrix, matrix), origin_matrix)


#
#  Other transformation routines
#


def transformation_matrix(
    translation: np.ndarray = np.zeros(3), rotation: np.ndarray = np.zeros(3)
) -> np.ndarray:
    """Returns a 4x4 transformation matrix."""
    matrix = rotation_matrix(rotation)
    matrix[:3, 3] = translation
    return matrix


def transform(coords: np.ndarray, matrix: np.ndarray):
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
    com = coords.mean(axis=0)

    vec1 = vector / np.linalg.norm(vector)
    vec2 = target / np.linalg.norm(target)

    axis = np.cross(vec1, vec2)
    sine = np.linalg.norm(axis)
    cosine = np.dot(vec1, vec2)
    amount = math.atan2(sine, cosine)

    return rotation_matrix_around_axis(axis, amount, center=com)


def orient(
    coords: np.ndarray,
    vector: np.ndarray | list[float],
    target: np.ndarray | list[float],
):
    """Orients coordinates.

    Args:
        vector (np.ndarray (3 x 1)): source vector
        target (np.ndarray (3 x 1)): target vector
    """
    t = orientation_matrix(coords, vector, target)
    transform(coords, t)


def dist(a: np.ndarray, b: np.ndarray) -> float:
    """Returns the euclidean distance between two sets of coordinates."""
    return (((a - b) ** 2.0).sum()) ** 0.5
