import os
from typing import Sequence, Union

import numpy as np


PathLike = Union[str, os.PathLike]
FilePath = PathLike
DirectoryPath = PathLike

Numeric = Union[float, int]
ArrayLike = Union[Sequence[Numeric], np.ndarray]