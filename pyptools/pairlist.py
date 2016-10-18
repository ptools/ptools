
"""pairlist - Iterate over atoms in interaction."""


import numpy
from scipy.spatial.distance import cdist


class PairList(list):
    """Iterate over atoms in interactions in two RigidBody instances.

    By construction, contacts stored in PairList are sorted by increasing
    atom indices.
    This sort applies for the second element of the pair as well, meaning that
    if a pair list contains [[1, 2], [1, 3]], element [1, 2] will always come
    before element [1, 3].

    Attrs:
        receptor (RigidBody)
        ligand (RigidBody)
        sqcutoff (float): squared cutoff for neighbor search

    Args:
        receptor (RigidBody)
        ligand (RigidBody)
        cutoff (float): cut-off for neighbor searching
    """
    def __init__(self, receptor, ligand, cutoff):
        super(PairList, self).__init__()
        self.receptor = receptor
        self.ligand = ligand
        self.sqcutoff = cutoff * cutoff
        self._contacts = []
        self.update()

    def update(self):
        """Update contact lists with the neighbor searching algorithm."""
        del self[:]
        dist = cdist(self.receptor.coords, self.ligand.coords, metric='sqeuclidean')
        idx = numpy.where(dist <= self.sqcutoff)
        self.extend(zip(*idx))

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
