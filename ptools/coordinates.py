from __future__ import annotations
from typing import Optional

import numpy as np
from numpy.typing import ArrayLike

import ptools.linalg as L


class Invalid3DCoordinates(Exception):
    """Exception raised when trying to initialize Coordinates3D from an invalid input array.

    Attributes:
        shape: tuple[int, ...] -- input array dimensions
    """

    def __init__(self, shape: tuple[int, ...]):
        message = f"cannot initialize 3D-coordinates from array with shape {shape}"
        super().__init__(message)


class Coordinates3D(np.ndarray):
    """Holds 3D-coordinates.

    Main purpose is to ensure the 3D-ness of the coordinates.
    """

    def __new__(cls, input_array: ArrayLike):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array, dtype="float64").view(cls)

        ndim = len(obj.shape)

        if not 1 <= ndim <= 2:
            raise Invalid3DCoordinates(obj.shape)
        if ndim == 1 and obj.shape[0] != 3:
            raise Invalid3DCoordinates(obj.shape)
        elif ndim == 2 and obj.shape[1] != 3:
            raise Invalid3DCoordinates(obj.shape)

        return obj

    def centroid(self) -> np.ndarray:
        """Returns an spatial object geometric center."""
        return L.centroid(self)

    def center(self, target: Optional[ArrayLike] = None):
        """Centers coordinates on target."""
        self = self.centered(target)

    def centered(self, target: Optional[ArrayLike] = None) -> Coordinates3D:
        """Returns the centered coordinates relatively to target."""
        if target is None:
            target = self.centroid()
        return self.__class__(self - target)

    def norm(self) -> float:
        return L.norm(self)

    def normalize(self):
        """Normalize coordinates."""
        self = self.normalized()

    def normalized(self) -> Coordinates3D:
        """Returns normalized coordinates."""
        return self.__class__(L.normalized(self))
