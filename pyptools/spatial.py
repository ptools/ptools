
"""pyptools.spatial - Defines classes and functions to work on spatial
(mostly 3D) data."""


import numpy


class SpatialObject(object):
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
        """Translate object coordinates using vector `v`."""
        self.coords += v


def coord3d(value=(0, 0, 0)):
    """Convert an object to a 1 x 3 shaped numpy array of floats."""
    if isinstance(value, (int, float)):
        return numpy.full((3,), value, dtype=float)
    value = numpy.array(value, dtype=float)
    if value.shape != (3, ):
        err = '3-d coordinates should be a scalar or '\
              '1 x 3 shaped-array ({} shaped-array found)'
        err = err.format(value.shape)
        raise ValueError(err)
    return value

