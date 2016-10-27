
"""pyptools.spatial - Defines classes and functions to work on spatial
(mostly 3D) data."""


import math

import numpy


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

    def translate(self, v):
        """Translate object coordinates using vector `v`.

        Args:
            v (array[float, int], int): 1 x 3 shaped vector or scalar
        """
        translate(self.coords, v)

    def rotate(self, alpha=0.0, beta=0.0, gamma=0.0):
        """Rotate object coordinates around X, Y and Z axes.

        Args:
            alpha (float): rotation angle around the X-axis
            beta (float): rotation angle around the Y-axis
            gamma (float): rotation angle around the Z-axis
        """
        rotate(self.coords, alpha, beta, gamma)

    def transform(self, matrix):
        """Transform object by 4x4 matrix.

        Args:
            matrix (numpy.ndarray): 4 x 4 matrix.
        """
        transform(self.coords, matrix)

    def attract_euler_rotate(self, phi, ssi, rot):
        """Rotate object with Attract convention."""
        attract_euler_rotate(self.coords, phi, ssi, rot)


def coord3d(value=(0, 0, 0)):
    """Convert an object to a 1 x 3 shaped numpy array of floats."""
    if isinstance(value, (int, float)):
        return numpy.full((3,), value, dtype=float)
    value = numpy.array(value, dtype=float)

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


def translate(coords, v):
    """In-place translation of coordinates by a vector.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        v (iterable[float, int]): 1 x 3 shaped vector
    """
    numpy.add(coords, v, coords)


def rotation(alpha=0.0, beta=0.0, gamma=0.0):
    """Return the rotation matrix around the X, Y and Z axes.

    The matrix rotates first along X-axis, then Y, then Z.

    Args:
        alpha (float): rotation angle around the X-axis (in degrees)
        beta (float): rotation angle around the Y-axis (in degrees)
        gamma (float): rotation angle around the Z-axis (in degrees)

    Returns:
        numpy.ndarray: 4x4 rotation matrix.
    """
    alpha = math.radians(alpha)
    beta = math.radians(beta)
    gamma = math.radians(gamma)
    r = numpy.identity(4)
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


def rotate(coords, alpha=0.0, beta=0.0, gamma=0.0):
    """In-place rotation of coordinates around X, Y and Z axes.

    Args:
        coords (numpy.ndarray): N x 3 shaped array
        alpha (float): rotation angle around the X-axis
        beta (float): rotation angle around the Y-axis
        gamma (float): rotation angle around the Z-axis
    """
    matrix = rotation(alpha, beta, gamma)
    coords[:] = numpy.inner(coords, matrix[:3, :3])


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
    a = numpy.ones((coords.shape[0], n + 1))
    a[:, :-1] = coords

    # Apply transformation matrix to coordinates.
    a[:] = numpy.inner(a, matrix)

    # Weight coordinates and reshape into N x 3 again.
    coords[:] = (a[:, 0:n].T / a[:, n]).T


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

    eulermat = numpy.identity(3)

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
    coords[:] = numpy.inner(coords, matrix)
