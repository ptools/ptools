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
from . import measure
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
        """Checks for equality."""
        if not isinstance(other, self.__class__):
            err = f"cannot compare {self.__class__.__qualname__} with object of type {type(other)}"
            raise TypeError(err)

        attrs = self.__class__.__dataclass_fields__.keys() - ("coordinates",)
        for name in attrs:
            if getattr(self, name) != getattr(other, name):
                return False

        return np.allclose(self.coordinates, other.coordinates)

    def copy(self) -> ObjectWithCoordinates:
        """Returns a copy of itself."""
        return self.__class__(self.coordinates.copy())

    def normalize(self):
        """Normalize coordinates."""
        self.coordinates.normalize()


class SupportsTranslation(ObjectWithCoordinates):
    """Object which coordinates can be translated."""

    def center_to_origin(self, origin: ArrayLike = np.zeros(3)):
        """Centers object on `origin`."""
        self.translate(np.asarray(origin) - measure.centroid(self))

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
        T.rotate_by(self.coordinates, angles)

    def rotate(self, rotation: ArrayLike):
        """Rotates object using rotation matrix or 3 angles.

        Args:
            rotation (np.ndarray): 3 x 3 matrix or 3 x 1 vector of angles.
        """
        T.rotate(self.coordinates, rotation)

    def ab_rotate(
        self, A: ArrayLike, B: ArrayLike, amount: float, degrees: bool = True
    ):
        """Rotates object using PTools rotation around axis."""
        T.ab_rotate(self.coordinates, A, B, amount, degrees)

    def attract_euler_rotate(self, angles: np.ndarray = np.zeros(3)):
        """Rotates object with Attract convention."""
        T.attract_euler_rotate(self.coordinates, angles)

    def orient(self, vector: ArrayLike, target: ArrayLike):
        """Orients a SpatialObject."""
        T.orient(self.coordinates, vector, target)


class SupportsTransformation(SupportsTranslation, SupportsRotation):
    """Object that can be translated and rotated."""

    def transform(self, matrix: ArrayLike):
        """Transforms object using 4x4 matrix."""
        T.transform(self.coordinates, matrix)

    def move(self, matrix: ArrayLike):
        """Transforms object using 4x4 matrix.

        This is an alias for `SpatialObject.transform`.
        """
        self.transform(matrix)
