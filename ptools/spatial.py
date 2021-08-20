"""ptools.spatial - Defines classes and functions to work on spatial
(mostly 3D) data."""


import math

import numpy as np
from scipy.linalg import expm, norm


class SpatialObject:
    """An object with coordinates.

    Implements basic spatial operations such as translation, rotation, etc.
    """

    def __init__(self, coords=(0, 0, 0)):
        self._coords = coord3d(coords)

    def dist(self, other):
        """Returns the euclidean distance between two objects."""
        return dist(self, other)

    def copy(self):
        """Returns a copy of itself."""
        return self.__class__(self.coords)

    @property
    def coords(self):
        """Get SpatialObject cartesian coordinates."""
        return self._coords

    @coords.setter
    def coords(self, pos):
        """Set SpatialObject cartesian coordinates."""
        self._coords = coord3d(pos)

    def centroid(self):
        """Returns an spatial object geometric center."""
        return centroid(self._coords)

    def translate(self, t):
        """Translate object coordinates using vector `t`.

        Args:
            t ((int, float) or np.ndarray): scalar, 1 x 3 shaped vector or
                4 x 4 matrix
        """
        translate(self.coords, t)

    def rotate_by(self, angles):
        """Rotate object coordinates around X, Y and Z axes.

        Args:
            angles ([float, float, float]): rotation angle around
                the X-, Y- and Z-axes
        """
        rotate_by(self.coords, angles)

    def rotate(self, rotation):
        """Rotate object using rotation matrix or 3 angles.

        Args:
            rotation (np.ndarray): 3 x 3 matrix of 3 x 1 vector of angles.
        """
        rotate(self.coords, rotation)

    def ab_rotate(self, A, B, amount):
        """PTools rotation around axis.

        Args:
            A (np.ndarray (3 x 1)): point A coordinates
            B (np.ndarray (3 x 1)): point B coordinates
            amount (float): angle of rotation
        """
        ab_rotate(self.coords, A, B, amount)

    def transform(self, matrix):
        """Transform object by 4x4 matrix.

        Args:
            matrix (numpy.ndarray): 4 x 4 matrix.
        """
        transform(self.coords, matrix)

    def move(self, matrix):
        """Transform object by 4x4 matrix.

        This is an alias for `SpatialObject.transform`.

        Args:
            matrix (numpy.ndarray): 4 x 4 matrix.
        """
        self.transform(matrix)

    def moveby(self, direction):
        """Translate object coordinates using vector `v`.

        This is an alias for `SpatialObject.translate`.

        Args:
            v (array[float, int], int): 1 x 3 shaped vector or scalar
        """
        self.translate(direction)

    def attract_euler_rotate(self, phi, ssi, rot):
        """Rotate object with Attract convention."""
        attract_euler_rotate(self.coords, phi, ssi, rot)

    def orient(self, vector, target):
        """Orients an AtomCollection.

        Args:
            vector (np.ndarray (3 x 1)): source vector
            target (np.ndarray (3 x 1)): target vector
        """
        orient(self.coords, vector, target)

    def tensor_of_inertia(self, sort=False):
        """Returns an AtomCollection inertial tensors.

        Args:
            sort (bool): if True, tensors are sorted by amplitude.
        """
        return tensor_of_inertia(self.coords, sort, method="fast")

    def distance_to_axis(self, axis):
        """Returns the SpatialObject distance to an arbitrary axis.

        Args:
            axis (np.array(3)): axis in 3-D
        """
        return np.linalg.norm(np.cross(self.coords, axis))


def coord3d(value=(0, 0, 0), *args):
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


def _coordinates_from_scalar(value):
    """Returns a numpy array of size 3 filled with value."""
    return np.full((3,), value, dtype=float)


def _coordinates_from_vector(value):
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


def centroid(x):
    """Return x centroid.

    Centroid is the average coordinate along axis 0.
    """
    return np.mean(x, axis=0)


def center_of_mass(coords, weights):
    """Return the center of mass (barycenter)."""
    weights = weights.reshape(-1, 1)
    return (coords * weights).sum(axis=0) / weights.sum()


def angle(u, v):
    """Returns the angle between two vectors in radians."""
    return math.acos(np.dot(u, v) / (norm(u) * norm(v)))


def _tensor_of_inertia_accurate(coords, weights):
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


def _tensor_of_inertia_fast(coords):
    """Returns the inertia tensors of a set of atoms.

    Fast: does not take into account the weights
    """
    com = centroid(coords)
    X = coords - com
    tensors = np.dot(X.transpose(), X)
    return tensors


def tensor_of_inertia(coords, weights=None, method="accurate"):
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


def principal_axes(tensor, sort=True):
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
def translation_matrix(direction=np.zeros(3)):
    """Return the matrix to translate by direction vector."""
    matrix = np.identity(4)
    matrix[:3, 3] = direction[:3]
    return matrix


def translate(coords, t):
    """In-place translation of coordinates by a vector.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        t ((int, float) or np.ndarray): scalar, 1 x 3 shaped vector or
            4 x 4 matrix
    """

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
def rotation_matrix(angles=np.zeros(3)):
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


def rotate_by(coords, angles=np.zeros(3)):
    """In-place rotation of coordinates around X, Y and Z axes.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        angles ([float, float, float]): rotation angle around
                the X-, Y- and Z-axes
    """
    matrix = rotation_matrix(angles)
    rotate(coords, matrix)


def rotate(coords, r):
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


def attract_euler_rotation_matrix(phi, ssi, rot):
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


def attract_euler_rotate(coords, phi, ssi, rot):
    """In-place Euler rotation of coordinates with the Attract convention.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        phi (float):
        ssi (float):
        rot (float):
    """
    matrix = attract_euler_rotation_matrix(phi, ssi, rot)
    coords[:] = np.inner(coords, matrix)


def ab_rotation_matrix(A, B, amount):
    """Returns the rotation matrix to rotate around axis (A, B) by amount (in radians)."""
    return rotation_matrix_around_axis(B - A, amount, A)


def ab_rotate(coords, A, B, amount):
    """Rotates coords around axis (A, B) by amount theta (in radians)."""
    matrix = rotation_matrix_around_axis(B - A, amount, A)
    transform(coords, matrix)


def rotation_matrix_around_axis(axis, amount, center=np.zeros(3)):
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


def transformation_matrix(translation=np.zeros(3), rotation=np.zeros(3)):
    """Returns a transformation matrix.

    Args:
        translation (np.ndarray(3)): translation components in x,y,z
        rotation (np.ndarray(3)): rotation components along x,y,z

    Returns:
        np.ndarray(4, 4)
    """
    matrix = rotation_matrix(rotation)
    matrix[:3, 3] = translation
    return matrix


def transform(coords, matrix):
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


def orientation_matrix(coords, vector, target):
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


def orient(coords, vector, target):
    """Orients coordinates.

    Args:
        vector (np.ndarray (3 x 1)): source vector
        target (np.ndarray (3 x 1)): target vector
    """
    t = orientation_matrix(coords, vector, target)
    transform(coords, t)


def dist(a, b):
    """Returns the distance between two SpatialObject instances."""
    return (((a.coords - b.coords) ** 2.0).sum()) ** 0.5
