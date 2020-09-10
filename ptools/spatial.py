
"""ptools.spatial - Defines classes and functions to work on spatial
(mostly 3D) data."""


import math

import numpy as np
from scipy.linalg import expm, norm
from scipy.spatial.transform import Rotation


class SpatialObject:
    """A object which coordinates.

    Implements basic spatial operations such as translation, rotation, etc.
    """
    def __init__(self, coords=(0, 0, 0)):
        self._coords = coord3d(coords)

    @property
    def coords(self):
        """Get atom cartesian coordinates."""
        return self._coords

    @coords.setter
    def coords(self, pos):
        """Set atom cartesian coordinates."""
        self._coords = coord3d(pos)

    def centroid(self):
        return centroid(self._coords)

    def translate(self, t):
        """Translate object coordinates using vector `t`.

        Args:
            t ((int, float) or np.array): scalar, 1 x 3 shaped vector or
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

    def rotate(self, r):
        """Rotate object using rotation matrix or 3 angles."""
        rotate(self.coords, r)

    def ab_rotate(self, A, B, angle):
        """PTools rotation around axis."""
        ab_rotate(self.coords, A, B, angle)

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


def coord3d(value=(0, 0, 0), *args):
    """Convert an object to a 1 x 3 shaped numpy array of floats."""
    def bad_coord3d_initialization():
        args_s = ", ".join(f"{str(a)}" for a in args)
        raise ValueError(f"Bad coord3d initialization: coord3d({value}, {args_s})")

    if args:
        if not isinstance(value, (int, float)):
            bad_coord3d_initialization()
        elif len(args) != 2:
            bad_coord3d_initialization()
        value = (value, args[0], args[1])

    if isinstance(value, (int, float)):
        return np.full((3,), value, dtype=float)
    value = np.array(value, dtype=float)

    if len(value.shape) == 1:
        if value.shape[0] != 3:
            err = '3-d coordinates should be a scalar or '\
                  '1 x 3 shaped-array ({} shaped-array found)'
            err = err.format(value.shape)
            raise ValueError(err)
    elif len(value.shape) == 2:
        if value.shape[1] != 3:
            err = '3-d coordinate array should be N x 3 '\
                  '({} found)'.format(value.shape)
            raise ValueError(err)
    else:
        err = '3-d coordinate array should have at most 2 dimensions '\
              '({} found)'.format(len(value.shape))
        raise ValueError(err)
    return value


def centroid(x):
    """Return x centroid.

    Centroid is the average coordinate along axis 0.
    """
    return np.mean(x, axis=0)


def norm(u):
    """Returns the norm of a vector."""
    return np.linalg.norm(u)


def angle(u, v):
    """Returns the angle between two vectors in radians."""
    return math.acos(np.dot(u, v) / (norm(u) * norm(v)))


#
# Translation routines
#
def translation_matrix(direction=[0, 0, 0]):
    """Return the matrix to translate by direction vector."""
    m = np.identity(4)
    m[:3, 3] = direction[:3]
    return m


def translate(coords, t):
    """In-place translation of coordinates by a vector.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        t ((int, float) or np.array): scalar, 1 x 3 shaped vector or
            4 x 4 matrix
    """
    def _translate_scalar(x):
        np.add(coords, t, coords)

    def isscalar(s):
        return isinstance(t, (float, int))

    if isscalar(t):
        _translate_scalar(t)
    else:
        t = np.array(t)
        if t.shape == (4, 4):
            t = t[:3, 3]
        elif t.shape != (3, ):
            raise ValueError("Dimensions error: expected 3 x 1 or 4 x 4 "
                            f"(got {t.shape})")
        np.add(coords, t, coords)



#
#  Rotation routines
#
def rotation_matrix(angles=[0.0, 0.0, 0.0]):
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

    r[1, 0] = math.cos(alpha) * math.sin(gamma) + math.sin(alpha) * math.sin(beta) * math.cos(gamma)
    r[1, 1] = math.cos(alpha) * math.cos(gamma) - math.sin(alpha) * math.sin(beta) * math.sin(gamma)
    r[1, 2] = -math.sin(alpha) * math.cos(beta)

    r[2, 0] = math.sin(alpha) * math.sin(gamma) - math.cos(alpha) * math.sin(beta) * math.cos(gamma)
    r[2, 1] = math.sin(alpha) * math.cos(gamma) + math.cos(alpha) * math.sin(beta) * math.sin(gamma)
    r[2, 2] = math.cos(alpha) * math.cos(beta)
    return r


def rotate_by(coords, angles=[0.0, 0.0, 0.0]):
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
        raise ValueError("Dimensions error: expected 3 x 1 or 4 x 4 "
                        f"(got {r.shape}) ")

    coords[:] = np.inner(coords, matrix[:3, :3])


# TODO: rename to attract_euler_rotation_matrix
def attract_euler_rotation(phi, ssi, rot):
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

    eulermat = np.identity(3)

    eulermat[0][0] = crot * cscp + srot * sp
    eulermat[0][1] = srot * cscp - crot * sp
    eulermat[0][2] = sscp

    eulermat[1][0] = crot * cssp - srot * cp
    eulermat[1][1] = srot * cssp + crot * cp
    eulermat[1][2] = sssp

    eulermat[2][0] = -crot * ss
    eulermat[2][1] = -srot * ss
    eulermat[2][2] = cs

    return eulermat


def attract_euler_rotate(coords, phi, ssi, rot):
    """In-place Euler rotation of coordinates with the Attract convention.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        phi (float):
        ssi (float):
        rot (float):
    """
    matrix = attract_euler_rotation(phi, ssi, rot)
    coords[:] = np.inner(coords, matrix)


def ab_rotation_matrix(A, B, angle):
    """Returns the rotation matrix to rotate around axis (A, B) by angle (in radians)."""
    return rotation_matrix_around_axis(B - A, angle)


def ab_rotate(coords, A, B, angle):
    """Rotates coords around axis (A, B) by angle theta (in radians)."""
    m = rotation_matrix_around_axis(B - A, angle)
    rotate(coords, m)


def rotation_matrix_around_axis(axis, angle):
    """Returns the rotation matrix to rotate around the axis by given angle.

    Args:
        axis (np.array): N x 3
        angle (float): angle in radians

    Returns:
        numpy.array: 4 x 4 rotation matrix
    """
    r = expm(np.cross(np.identity(3), axis / norm(axis) * angle))
    m = np.identity(4)
    m[:3, :3] = r
    return m




#
#  Other transformation routines
#

def transformation_matrix(translation=[0., 0., 0.], rotation=[0., 0., 0.]):
    """Returns a transformation matrix.

    Args:
        translation (np.array(3)): translation components in x,y,z
        rotation (np.array(3)): rotation components along x,y,z

    Returns:
        np.array(4, 4)
    """
    m = rotation_matrix(rotation)
    m[:3, 3] = translation
    return m


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


