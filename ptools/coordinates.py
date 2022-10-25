from __future__ import annotations
from dataclasses import dataclass
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



class Coordinates3D:
    """Holds 3D-coordinates.

    More or less compatible with numpy.ndarray.
    Main purpose is to ensure the 3D-ness of the coordinates.
    """

    array: np.ndarray = np.zeros((0, 3))

    def __init__(self, array: ArrayLike = np.zeros((0, 3))):
        self.array = self._to_3d_array(array)

    def __repr__(self) -> str:
        class_name = f"{self.__class__.__qualname__}"
        if self.shape == (3, ):
            return f"{class_name}({self.array.tolist()!r})"

        if self.shape[0] < 5:
            repr_array = ", ".join(repr(row.tolist()) for row in self.array)
            return f"{class_name}([{repr_array}])"

        sep = f"\n{' ' * (len(class_name) + 2)}"
        rep = (
            f"{class_name}(["
            f"{sep.join(repr(row.tolist()) for row in self.array[:2])}"
            f"{sep}..."
            f"{sep}"
            f"{sep.join(repr(row.tolist()) for row in self.array[-2:])}"
            "])"
        )
        return rep

    @staticmethod
    def _to_3d_array(array: ArrayLike) -> np.ndarray:
        """Returns a valid input array, i.e. representing 3D-coordinates.

        Input array should be either in form (x y z) for ((x1 y1 z1), ..., (xn yn zn)).

        Raises:
            Invalid3DCoordinates: if array as incompatible dimensions
        """
        array = np.asarray(array, dtype=float)
        ndim = len(array.shape)

        if not 1 <= ndim <= 2:
            raise Invalid3DCoordinates(array.shape)
        if ndim == 1 and array.shape[0] != 3:
            raise Invalid3DCoordinates(array.shape)
        elif ndim == 2 and array.shape[1] != 3:
            raise Invalid3DCoordinates(array.shape)

        return array

    @property
    def shape(self) -> tuple[int, ...]:
        return self.array.shape

    def __eq__(self, __o: object) -> bool:
        if isinstance(__o, type(self)):
            return np.allclose(self.array, __o.array)
        return np.allclose(np.asarray(__o), self.array)

    def __getitem__(self, i: slice|int) -> Coordinates3D:
        return self.__class__(self.array[i])

    def __setitem__(self, i: slice|int, item):
        self.array[i] = item

    def __array__(self, *args) -> np.ndarray:
        return np.asarray(self.array, *args)

    def centroid(self) -> np.ndarray:
        """Returns an spatial object geometric center."""
        return L.centroid(self.array)

    def center(self, target: Optional[ArrayLike] = None):
        """Centers coordinates on target."""
        self = self.centered(target)

    def centered(self, target: Optional[ArrayLike] = None) -> Coordinates3D:
        """Returns the centered coordinates relatively to target."""
        if target is None:
            target = self.centroid()
        return self.__class__(self.array - target)

    def norm(self) -> float:
        return L.norm(self.array)

    def normalize(self):
        """Normalize coordinates."""
        self.array = L.normalized(self.array)

    def normalized(self) -> Coordinates3D:
        """Returns normalized coordinates."""
        return self.__class__(L.normalized(self.array))

