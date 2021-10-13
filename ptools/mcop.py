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
