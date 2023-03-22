import numpy as np
from numpy.typing import DTypeLike

from ptools.particlecollection import ParticleCollection
from .measure import bounding_box

class Grid:

    coordinates: np.ndarray
    data: np.ndarray
    _dimensions: tuple[int, int, int]

    def __init__(self, grid: np.ndarray, dtype: DTypeLike):
        self.coordinates = grid
        self.spacing = np.linalg.norm(grid[1] - grid[0])
        self.data = np.zeros(grid.shape[0], dtype=dtype)
        self._dimensions = (0, 0, 0)

    def __repr__(self) -> str:
        return f"{self.__class__.__qualname__}({self.coordinates}, {self.data.dtype})"

    @classmethod
    def from_boundaries(cls, bounds: np.ndarray, spacing: float, dtype: DTypeLike):
        """Generate grid"""
        x = np.arange(bounds[0, 0], bounds[1, 0] + spacing, spacing)
        y = np.arange(bounds[0, 1], bounds[1, 1] + spacing, spacing)
        z = np.arange(bounds[0, 2], bounds[1, 2] + spacing, spacing)
        grid = cls(np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3), dtype)
        grid._dimensions = (len(x), len(y), len(z))
        return grid

    @property
    def shape(self) -> tuple[int, int, int]:
        """Return grid shape"""
        return self._dimensions

    @property
    def origin(self) -> np.ndarray:
        """Returns the grid's origin"""
        return self.coordinates[0]

    def size(self) -> int:
        """Returns the number of points in the grid.

        Equivalent to ``np.prod(self.shape)``.
        """
        return self.coordinates.shape[0]

    def __iter__(self):
        """Iterate over grid points and data."""
        return iter(zip(self.coordinates, self.data))



def generate_grid(bounds: np.ndarray, spacing: float, dtype: DTypeLike) -> Grid:
    """Generate grid"""
    return Grid.from_boundaries(bounds, spacing, dtype)


def generate_solvation_grid(atoms: ParticleCollection, spacing: float) -> Grid:
    """Generates a grid where points that are not within the radius of any atom are marked as solvent."""
    grid = generate_grid(bounding_box(atoms), spacing, bool)
    grid.data = np.ones(grid.size(), dtype=bool)
    for atom in atoms:
        distance = np.linalg.norm(grid.coordinates - atom.coordinates, axis=1)
        grid.data[(distance - atom.radius) < 1e-6] = False
    return grid

