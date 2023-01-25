"""test_pairlist - Tests for `ptools.pairlist` module."""

import unittest

import numpy as np

import ptools
from ptools.pairlist import PairList

from . import TEST_LIGAND, TEST_RECEPTOR, TEST_DISTANCES_RECEPTOR_LIGAND
from .testing.moreassert import assert_array_almost_equal




class TestPairList(unittest.TestCase):
    def setup_class(self):
        self.cutoff = 5.0
        self.ligand = ptools.read_pdb(TEST_LIGAND)
        self.receptor = ptools.read_pdb(TEST_RECEPTOR)
        assert len(self.ligand) == 974
        assert len(self.receptor) == 936
        self.pairlist = PairList(self.receptor, self.ligand, self.cutoff)

    def test_contacts(self):
        """Checks that contacts measured with PairList are the same as the ones
        measured with VMD."""
        actual = self.pairlist.contacts()
        expected = self.read_reference_contacts()
        assert actual == expected

    def test_distances(self):
        """Checks that the number of distances corresponds to the number of
        contacts and that all distances are below cutoff."""
        distances = self.pairlist.distances()
        assert len(distances) == len(self.pairlist.contacts())
        assert all(distances < self.cutoff)

    def test_all_distances(self):
        expected = self.read_reference_distance_file()
        distances = self.pairlist.all_distances()
        assert distances.shape == expected.shape

        # approx doesn't support nested data structures
        assert_array_almost_equal(distances, expected, decimal=5)

    def read_reference_distance_file(self) -> np.ndarray:
        """Reads the file that contains all the distances between `TEST_RECEPTOR` and `TEST_LIGAND`."""
        data = np.loadtxt(TEST_DISTANCES_RECEPTOR_LIGAND)
        assert data.shape == (self.receptor.size(), self.ligand.size())
        return data

    def read_reference_contacts(self) -> list[tuple[int, int]]:
        distances = self.read_reference_distance_file()
        indices = np.where(distances < self.cutoff)
        return list(zip(*indices))

