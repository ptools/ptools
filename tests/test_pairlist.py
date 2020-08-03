
"""test_pairlist - Tests for `ptools.pairlist` module."""

import unittest

from ptools.rigidbody import RigidBody
from ptools.pairlist import PairList

from . import TEST_LIGAND, TEST_RECEPTOR, TEST_DISTANCES_RECEPTOR_LIGAND
from .testing.moreassert import assert_array_equal, assert_array_almost_equal


# Contacts between receptor.pdb and ligand.pdb from tests/data.
# Contacts have been measured with VMD with a cut-off of 5.0 Å.
# CONTACTS_REF is two lists of atom indices, the first containing the first
# index of each pair and the second containing the second index.
CONTACTS_REF = [
    [843, 842, 843, 842, 834, 343, 824, 835, 834, 848, 850, 851, 851, 847,
     849, 856, 337, 342, 343, 345, 341, 343, 851, 436, 435, 410, 341, 859,
     858, 859, 858, 436, 435, 410, 436, 443, 163, 164, 163, 164, 163, 165,
     164, 176, 160, 442, 443, 22, 160, 160, 30, 30, 97, 153, 154],
    [973, 773, 773, 774, 775, 781, 778, 778, 776, 567, 567, 567, 568, 567,
     567, 572, 780, 561, 561, 561, 561, 560, 560, 734, 724, 724, 786, 558,
     558, 559, 559, 727, 723, 723, 723, 721, 551, 551, 552, 552, 553, 552,
     548, 550, 550, 720, 720, 541, 288, 287, 536, 537, 528, 283, 283]
]


class TestPairList(unittest.TestCase):

    def setUp(self):
        self.cutoff = 5.0
        self.ligand = RigidBody(TEST_LIGAND)
        self.receptor = RigidBody(TEST_RECEPTOR)
        self.assertEqual(len(self.ligand), 974)
        self.assertEqual(len(self.receptor), 936)
        self.pairlist = PairList(self.receptor, self.ligand, self.cutoff)

    def test_contacts(self):
        # Check that contacts measured with PairList are the same as
        # those the reference (calculated from VMD).
        pairs = self.pairlist.contacts()
        pairs_ref = PairList.sort(*CONTACTS_REF)
        assert_array_equal(pairs, pairs_ref)

    def test_distances(self):
        # Check that the number of distances corresponds to the number of
        # contacts and that all distances are below cutoff.
        distances = self.pairlist.distances()
        self.assertEqual(len(distances), len(self.pairlist.contacts()))
        self.assertTrue(all(distances < self.cutoff))

    def test_all_distances(self):
        # Read reference distances.
        with open(TEST_DISTANCES_RECEPTOR_LIGAND, 'rt') as f:
            ref = [[float(value) for value in line.split()]
                   for line in f if not line.startswith('#')]
        distances = self.pairlist.all_distances()
        self.assertEqual(distances.shape, (len(ref), len(ref[0])))
        assert_array_almost_equal(distances, ref, decimal=5)
