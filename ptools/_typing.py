import os
from collections.abc import Sequence
from typing import Protocol

import numpy as np

PathLike = str | os.PathLike
FilePath = PathLike
DirectoryPath = PathLike

Numeric = float | int
ArrayLike = Sequence[Numeric] | np.ndarray


class HasCoordinatesType(Protocol):
    @property
    def coordinates(self) -> np.ndarray:
        ...


class TopologyType(Protocol):
    @property
    def coordinates(self) -> np.ndarray:
        ...

    @property
    def masses(self) -> np.ndarray:
        ...


class AtomType(Protocol):
    @property
    def name(self) -> str:
        ...

    @property
    def index(self) -> int:
        ...

    @property
    def residue_name(self) -> str:
        ...

    @property
    def residue_index(self) -> int:
        ...

    @property
    def chain(self) -> str:
        ...

    @property
    def coordinates(self) -> ArrayLike:
        ...
