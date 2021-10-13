"""Defines classes and functions used for the multicopy feature."""

from typing import Sequence
from dataclasses import dataclass, field

from .atom import AtomCollection


@dataclass
class Mcop:
    """Container for multiple copies of the same monomer."""

    copies: Sequence[AtomCollection] = field(default_factory=list)

    def add_copy(self, kopy: AtomCollection):
        """Append a copy to the list of copies."""
        self.copies.append(kopy)

    def size(self):
        """Returns the number of copies.

        Alias:
            len(mcop)
        """
        return len(self.copies)

    def __len__(self):
        """Returns the number of copies.

        Alias:
            mcop.size()
        """
        return self.size()

    def __getitem__(self, key) -> AtomCollection:
        """Returns the copy at key `key`."""
        return self.copies[key]

    def clear(self):
        """Clears the internal list of copies."""
        self.copies.clear()