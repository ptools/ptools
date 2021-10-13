"""Tests for ptools.mcop"""

import unittest

from ptools.atom import AtomCollection, BaseAtom
from ptools.mcop import Mcop

from .testing.moreassert import assert_array_almost_equal


class TestMcop(unittest.TestCase):
    """Tests for the Mcop class."""

    def test_add_copy(self):
        mcop = Mcop()
        self.assertEqual(len(mcop.copies), 0)
        mcop.add_copy(dummy_atom_collection())
        self.assertEqual(len(mcop.copies), 1)

    def test_add_copy_makes_copies(self):
        col = dummy_atom_collection()
        mcop = Mcop()
        mcop.add_copy(col)
        col.coords += 12  # after modification, mcop copy should be different.
        self.assertFalse(assert_array_almost_equal(col.coords, mcop.copies[0].coords))



def dummy_atom_collection(n_atoms=10):
    """Create a dummy AtomCollection composed of `n_atoms` atoms."""
    return AtomCollection([BaseAtom(coords=(i, i, i)) for i in range(n_atoms)])
