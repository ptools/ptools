"""ptools.spatial - Defines classes and functions to work on spatial
(mostly 3D) data."""

# Allows to use type hinting of class within itself
# e.g. SpatialObject.dist(self, other: SpatialObject)
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .array3d import array3d
from ._typing import ArrayLike

from . import linalg as L
from .linalg import transform as T


def default_array3d():
    return array3d(np.zeros((3,)))


@dataclass
class ObjectWithCoordinates:
    """Object with coordinates.

    Coordinates are stored in private _coords array which allows to automatically
    convert coordinates to numpy arrays upon change.
    """

    coordinates: array3d = field(default_factory=default_array3d)

    def __post_init__(self):
        self.coordinates = array3d(self.coordinates)

    def __eq__(self, other: object) -> bool:
        """Compares two BaseAtom instances."""
        if not isinstance(other, self.__class__):
            err = f"cannot compare {self.__class__.__qualname__} with object of type {type(other)}"
            raise TypeError(err)

        attrs = self.__class__.__dataclass_fields__.keys() - ("coordinates",)
        for name in attrs:
            if getattr(self, name) != getattr(other, name):
                return False

        return np.allclose(self.coordinates, other.coordinates)

    @property
    def coords(self) -> array3d:
        """Returns SpatialObject cartesian coordinates."""
        return self.coordinates

    @coords.setter
    def coords(self, pos: ArrayLike):
        """Sets SpatialObject cartesian coordinates."""
        self.coordinates = array3d(pos)

    def copy(self) -> ObjectWithCoordinates:
        """Returns a copy of itself."""
        return self.__class__(self.coordinates.copy())

    def distance(self, other: ObjectWithCoordinates) -> float:
        """Returns the euclidean distance between two objects."""
        return L.distance(self.coords, other.coords)

    def centroid(self) -> np.ndarray:
        """Returns an spatial object geometric center."""
        return self.coordinates.centroid()

    def distance_to_axis(self, axis: ArrayLike) -> float:
        """Returns the SpatialObject distance to an arbitrary axis."""
        return L.distance_to_axis(self.coords, axis)

    def normalize(self):
        """Normalize coordinates."""
        self.coordinates.normalize()


class SupportsTranslation(ObjectWithCoordinates):
    """Object which coordinates can be translated."""

    def center_to_origin(self, origin: ArrayLike = np.zeros(3)):
        """Centers object on `origin`."""
        self.translate(np.asarray(origin) - self.centroid())

    def translate(self, direction: ArrayLike):
        """Translates object coordinates using vector `direction`.

        Args:
            direction (scalar or array): scalar, 1 x 3 shaped vector or
                4 x 4 matrix
        """
        T.translate(self.coordinates, direction)

    def moveby(self, direction: ArrayLike):
        """Translates object coordinates using vector or scalar `direction`.

        This is an alias for `SpatialObject.translate`.
        """
        self.translate(direction)


class SupportsRotation(ObjectWithCoordinates):
    """Object which coordinates can be rotated."""

    def rotate_by(self, angles: ArrayLike):
        """Rotates object coordinates around X, Y and Z axes.

        Args:
            angles (3, ): rotation angles around the X-, Y- and Z-axes
        """
        T.rotate_by(self.coords, angles)

    def rotate(self, rotation: ArrayLike):
        """Rotates object using rotation matrix or 3 angles.

        Args:
            rotation (np.ndarray): 3 x 3 matrix or 3 x 1 vector of angles.
        """
        T.rotate(self.coords, rotation)

    def ab_rotate(
        self, A: ArrayLike, B: ArrayLike, amount: float, degrees: bool = True
    ):
        """Rotates object using PTools rotation around axis."""
        T.ab_rotate(self.coords, A, B, amount, degrees)

    def attract_euler_rotate(self, angles: np.ndarray = np.zeros(3)):
        """Rotates object with Attract convention."""
        T.attract_euler_rotate(self.coords, angles)

    def orient(self, vector: ArrayLike, target: ArrayLike):
        """Orients a SpatialObject."""
        T.orient(self.coords, vector, target)


class SupportsTransformation(SupportsTranslation, SupportsRotation):
    """Object that can be translated and rotated."""

    def transform(self, matrix: ArrayLike):
        """Transforms object using 4x4 matrix."""
        T.transform(self.coords, matrix)

    def move(self, matrix: ArrayLike):
        """Transforms object using 4x4 matrix.

        This is an alias for `SpatialObject.transform`.
        """
        self.transform(matrix)


# pylint: disable-msg=W1113
# keyword-arg-before-vararg
def coord3d(value: ArrayLike = np.zeros(3), *args) -> np.ndarray:
    """Converts an iterable of size 3 to a 1 x 3 shaped numpy array of floats.

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
            raise ValueError(
                f"Coordinates must be initialized either "
                f"with 1 or 3 arguments (found {value=}, {args=})"
            )

    # Initialization from a single value.
    if not args:
        # value is a scalar.
        if np.isscalar(value):
            return _coordinates_from_scalar(_to_float(value))
        if isinstance(value, (np.ndarray, Sequence)):
            return _coordinates_from_vector(value)
        # value is some invalid type.
        raise TypeError(f"Invalid coordinates initialization from {value}")

    # Initialization from 3 values.
    return np.array((value, *args))


def _coordinates_from_scalar(value: int | float) -> np.ndarray:
    """Returns a numpy array of size 3 filled with value."""
    return np.full((3,), value, dtype=float)


def _coordinates_from_vector(value: ArrayLike) -> np.ndarray:
    """Returns a numpy array of size 3 initialized with a vector."""
    array = np.array(value, dtype=float)

    ndim = np.ndim(array)
    if ndim > 2:
        raise ValueError(
            "3D coordinate array should have at "
            f"most 2 dimensions (found {array.shape}"
        )
    if (ndim == 2 and array.shape[1] != 3) or (ndim == 1 and array.shape[0] != 3):
        raise ValueError(
            "3D coordinate array should be N x 3 " f"(found {array.shape})"
        )
    return array


def _to_float(arg: Any) -> float:
    try:
        return float(arg)
    except:
        raise
