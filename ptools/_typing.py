import os
from typing import Protocol, Sequence, Union

import numpy as np


PathLike = Union[str, os.PathLike]
FilePath = PathLike
DirectoryPath = PathLike

Numeric = Union[float, int]
ArrayLike = Union[Sequence[Numeric], np.ndarray]


class HasCoordinates(Protocol):
    @property
    def coordinates(self) -> ArrayLike:
        ...

class Topology(Protocol):
    @property
    def coordinates(self) -> ArrayLike:
        ...

    @property
    def masses(self) -> ArrayLike:
        ...
