from dataclasses import dataclass, field
from typing import Optional

import numpy as np


@dataclass
class Balloon:
    """Dummy object with coordinates."""

    coordinates: np.ndarray = field(default_factory=lambda: np.zeros((5, 3)))


def generate_balloon(coordinates: Optional[np.typing.ArrayLike] = None) -> Balloon:
    """Returns a ``Balloon`` with given coordinates"""
    if coordinates is None:
        return Balloon()
    return Balloon(np.asarray(coordinates))
