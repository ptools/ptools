
"""pyptools.spatial - Defines classes and functions to work on spatial
(mostly 3D) data."""


import math
from math import cos, sin

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
    r[0, 0] = cos(beta) * cos(gamma)
    r[0, 1] = -cos(beta) * sin(gamma)
    r[0, 2] = sin(beta)

    r[1, 0] = cos(alpha) * sin(gamma) + sin(alpha) * sin(beta) * cos(gamma)
    r[1, 1] = cos(alpha) * cos(gamma) - sin(alpha) * sin(beta) * sin(gamma)
    r[1, 2] = -sin(alpha) * cos(beta)

    r[2, 0] = sin(alpha) * sin(gamma) - cos(alpha) * sin(beta) * cos(gamma)
    r[2, 1] = sin(alpha) * cos(gamma) + cos(alpha) * sin(beta) * sin(gamma)
    r[2, 2] = cos(alpha) * cos(beta)
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
    return numpy.dot(coords, matrix[:3, :3])



