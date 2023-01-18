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
    "coordinates": (1, 2, 3),
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
    from ptools.atomattrs import guess_atom_element, guess_atom_mass

    atoms = [AtomAttrs(coordinates=(i, i, i)) for i in range(size)]
    for atom in atoms:
        for key, value in DUMMY_ATOM_ATTRS.items():
            if key != "coords":
                setattr(atom, key, value)
        atom.guess_element()
        atom.guess_mass()

    return AtomCollection.from_objects(atoms)
