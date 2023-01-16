"""Creates dummy ptools structures."""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ptools.atomattrs import AtomAttrs
from ptools.atomcollection import AtomCollection

DUMMY_ATOM_ATTRS = {
    "name": "CA",
    "index": 42,
    "residue_name": "ALA",
    "residue_index": 17,
    "chain": "A",
    "charge": 2.0,
    "coords": (1, 2, 3),
}


@dataclass
class Balloon:
    """Dummy object with coordinates."""

    coordinates: np.ndarray = field(default_factory=lambda: np.zeros((5, 3)))


def generate_balloon(coordinates: Optional[np.typing.ArrayLike] = None) -> Balloon:
    """Returns a ``Balloon`` with given coordinates"""
    if coordinates is None:
        return Balloon()
    return Balloon(np.asarray(coordinates))


def generate_dummy_atom() -> AtomAttrs:
    """Creates a dummy atom."""
    atom = AtomAttrs()
    for attr, value in DUMMY_ATOM_ATTRS.items():
        setattr(atom, attr, value)
    return atom


def generate_dummy_atomcollection(size: int = 10) -> AtomCollection:
    """Creates a dummy atom collection composed of `size` atoms.

    Atom coordinates are [(0, 0, 0), (1, 1, 1), ..., (9, 9, 9)].
    """
    col = AtomCollection([AtomAttrs(coordinates=(i, i, i)) for i in range(size)])
    for atom in col:
        for attr in (
            "name",
            "index",
            "residue_name",
            "residue_index",
            "chain",
            "charge",
        ):
            setattr(atom, attr, DUMMY_ATOM_ATTRS[attr])
    col.guess_masses()
    return col
