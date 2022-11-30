"""pairlist - Iterate over atoms in interaction."""

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial.distance import cdist

from .atomcollection import AtomCollection


def zeros3f():
    return np.zeros((3,), dtype="float64")


@dataclass
class PairList:
    """Retrieve atoms that are within a certain cutoff in AtomCollections.

    By construction, contacts stored in PairList are sorted by increasing
    atom indices.
    This sort applies for the second element of the pair as well, meaning that
    if a pair list contains [[1, 2], [1, 3]], element [1, 2] will always come
    before element [1, 3].
    """

    receptor: AtomCollection
    ligand: AtomCollection
    cutoff: float
    _all_sqdistances: np.ndarray = field(
        init=False, repr=False, default_factory=zeros3f
    )
    _contacts: np.ndarray = field(init=False, repr=False, default_factory=zeros3f)

    def __post_init__(self):
        self.update()

    def contacts(self) -> list[tuple[int, int]]:
        """Get the indices of atoms pairs that are within a cutoff of
        each other."""
        return list(zip(*self._contacts))

    def all_sqdistances(self) -> np.ndarray:
        """Return the matrix of squared distances between every atom pairs."""
        return self._all_sqdistances

    def sqdistances(self) -> np.ndarray:
        """Return the matrix of squared distances between atoms within cutoff."""
        return self._all_sqdistances[self._contacts]

    def all_distances(self) -> np.ndarray:
        """Return the matrix of distances between every atom pairs."""
        return np.sqrt(self.all_sqdistances())

    def distances(self) -> np.ndarray:
        """Return the matrix of distances between atoms within cutoff."""
        return np.sqrt(self.sqdistances())

    def update(self):
        """Update contact and distance lists with the neighbor searching
        algorithm."""
        self._all_sqdistances = cdist(
            self.receptor.coordinates, self.ligand.coordinates, metric="sqeuclidean"
        )
        sqcutoff = self.cutoff * self.cutoff
        self._contacts = np.where(self._all_sqdistances <= sqcutoff)

    @staticmethod
    def sort(receptor: list[int], ligand: list[int]) -> list[tuple[int, int]]:
        """Sort a contact list.

        Returns:
            list[(int, int)]: pairs of atom indexes.
        """
        return sorted(zip(receptor, ligand), key=lambda x: (x[0], x[1]))
