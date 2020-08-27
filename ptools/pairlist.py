
"""pairlist - Iterate over atoms in interaction."""


import numpy
from scipy.spatial.distance import cdist


class PairList:
    """Retrieve atoms that are within a certain cutoff in AtomCollections.

    By construction, contacts stored in PairList are sorted by increasing
    atom indices.
    This sort applies for the second element of the pair as well, meaning that
    if a pair list contains [[1, 2], [1, 3]], element [1, 2] will always come
    before element [1, 3].

    Attrs:
        receptor (AtomCollection)
        ligand (AtomCollection)
        sqcutoff (float): squared cutoff for neighbor search

    Args:
        receptor (AtomCollection)
        ligand (AtomCollection)
        cutoff (float): cut-off for neighbor searching
    """
    def __init__(self, receptor, ligand, cutoff):
        super(PairList, self).__init__()
        self.receptor = receptor
        self.ligand = ligand
        self.sqcutoff = cutoff * cutoff
        self._all_distances = None
        self._contacts = None
        self.update()

    def contacts(self):
        """Get the indices of atoms pairs that are within a cutoff of
        each other."""
        return list(zip(*self._contacts))

    def all_sqdistances(self):
        """Return the matrix of squared distances between every atom pairs."""
        return self._all_distances

    def sqdistances(self):
        """Return the matrix of squared distances between atoms within cutoff."""
        return self._all_distances[self._contacts]

    def all_distances(self):
        """Return the matrix of distances between every atom pairs."""
        return numpy.sqrt(self.all_sqdistances())

    def distances(self):
        """Return the matrix of distances between atoms within cutoff."""
        return numpy.sqrt(self.sqdistances())

    def update(self):
        """Update contact and distance lists with the neighbor searching
        algorithm."""
        self._all_distances = cdist(self.receptor.coords, self.ligand.coords,
                                    metric='sqeuclidean')
        self._contacts = numpy.where(self._all_distances <= self.sqcutoff)

    @staticmethod
    def sort(receptor, ligand):
        """Sort a contact list.

        Args:
            receptor (list[int]): receptor atom indexes
            ligand (list[int]): ligand atom indexes

        Returns:
            list[(int, int)]: pairs of atom indexes.
        """
        return sorted(zip(receptor, ligand), key=lambda x: (x[0], x[1]))