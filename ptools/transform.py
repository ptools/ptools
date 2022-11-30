"""ptools.transform - Moves objects around."""

import numpy as np

from . import measure
from .linalg import transform as T
from ._typing import ArrayLike, HasCoordinates


def translate(obj: HasCoordinates, direction: ArrayLike):
    """Translates an object coordinates using vector `direction`.

    Args:
        direction (scalar or array): scalar, 1 x 3 shaped vector or
            4 x 4 matrix
    """
    T.translate(obj.coordinates, direction)


def moveby(obj: HasCoordinates, direction: ArrayLike):
    """Translates object coordinates using vector or scalar `direction`.

    This is an alias for ``translate``.
    """
    translate(obj, direction)


def center_to_origin(obj: HasCoordinates, origin: ArrayLike = np.zeros(3)):
    """Centers object on `origin`."""
    translate(obj, np.asarray(origin) - measure.centroid(obj))


def rotate_by(obj: HasCoordinates, angles: ArrayLike):
    """Rotates object coordinates around X, Y and Z axes.

    Args:
        angles (3, ): rotation angles around the X-, Y- and Z-axes
    """
    T.rotate_by(obj.coordinates, angles)


def rotate(obj: HasCoordinates, rotation: ArrayLike):
    """Rotates object using rotation matrix or 3 angles.

    Args:
        rotation (np.ndarray): 3 x 3 matrix or 3 x 1 vector of angles.
    """
    T.rotate(obj.coordinates, rotation)


def ab_rotate(
    obj: HasCoordinates, A: ArrayLike, B: ArrayLike, amount: float, degrees: bool = True
):
    """Rotates object using PTools rotation around axis."""
    T.ab_rotate(obj.coordinates, A, B, amount, degrees)


def attract_euler_rotate(obj: HasCoordinates, angles: np.ndarray = np.zeros(3)):
    """Rotates object with Attract convention."""
    T.attract_euler_rotate(obj.coordinates, angles)


def orient(obj: HasCoordinates, vector: ArrayLike, target: ArrayLike):
    """Orients a SpatialObject."""
    T.orient(obj.coordinates, vector, target)


def transform(obj: HasCoordinates, matrix: ArrayLike):
    """Transforms an object using 4x4 matrix."""
    T.transform(obj.coordinates, matrix)


def move(obj: HasCoordinates, matrix: ArrayLike):
    """Transforms an object using 4x4 matrix.

    This is an alias for ``transform``.
    """
    transform(obj, matrix)
