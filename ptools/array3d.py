from __future__ import annotations
from typing import Any, Optional

import numpy as np
from numpy.typing import ArrayLike

import ptools.linalg as L


class Invalid3DArrayError(ValueError):
    """Exception raised when trying to initialize array3d from an invalid input array.

    Attributes:
        shape: tuple[int, ...] -- input array dimensions
    """

    def __init__(self, shape: tuple[int, ...]):
        message = f"cannot initialize 3D-coordinates from array with shape {shape}"
        super().__init__(message)


class array3d(np.ndarray):
    """Holds 3D-coordinates.

    Main purpose is to ensure the 3D-ness of the coordinates.
    """

    def __new__(cls, input_array: ArrayLike):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array, dtype="float64").view(cls)

        # Checks input_array dimensions.
        validate_array3d_dimensions(obj)

        return obj

    @classmethod
    def zeros(cls, shape: int | tuple[int, ...] = 3) -> array3d:
        """Builds a 0-filled 3D-array."""
        return cls(np.zeros(shape, dtype=float))

    def centroid(self) -> np.ndarray:
        """Returns an spatial object geometric center."""
        return L.centroid(self)

    def center(self, target: Optional[ArrayLike] = None):
        """Centers coordinates on target."""
        self = self.centered(target)

    def centered(self, target: Optional[ArrayLike] = None) -> array3d:
        """Returns the centered coordinates relatively to target."""
        if target is None:
            target = self.centroid()
        return self.__class__(self - target)

    def norm(self) -> float:
        return L.norm(self)

    def normalize(self):
        """Normalize coordinates."""
        self = self.normalized()

    def normalized(self) -> array3d:
        """Returns normalized coordinates."""
        return self.__class__(L.normalized(self))


def asarray3d(obj: ArrayLike) -> array3d:
    """Converts an object to an ``array3d``."""
    return array3d(obj)


def validate_array3d_dimensions(obj: Any):
    """Raises Invalid3DArrayError if ``obj`` dimensions do not correspond to 3D-data."""
    shape = np.shape(obj)
    ndim = len(shape)

    if not 1 <= ndim <= 2:
        raise Invalid3DArrayError(shape)
    if ndim == 1 and shape[0] != 3:
        raise Invalid3DArrayError(shape)
    elif ndim == 2 and shape[1] != 3:
        raise Invalid3DArrayError(shape)
