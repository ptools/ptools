"""Creates dummy ptools structures."""

import numpy as np

from ptools.atom import BaseAtom, AtomCollection


DUMMY_ATOM_ATTRS = {
    "name": "CA",
    "resname": "ALA",
    "chain": "A",
    "index": 42,
    "resid": 17,
    "charge": 2.0,
    "coords": (1, 2, 3),
}


def generate_dummy_atom() -> BaseAtom:
    """Creates a dummy atom."""
    atom = BaseAtom()
    for attr, value in DUMMY_ATOM_ATTRS.items():
        setattr(atom, attr, value)
    return atom


def generate_dummy_atomcollection(size: int = 10) -> AtomCollection:
    """Creates a dummy atom collection composed of `size` atoms.

    Atom coordinates are [(0, 0, 0), (1, 1, 1), ..., (9, 9, 9)].
    """
    col = AtomCollection([BaseAtom(coords=(i, i, i)) for i in range(size)])
    for atom in col:
        for attr in ('name', 'resname', 'chain', 'index', 'resid', 'charge'):
            setattr(atom, attr, DUMMY_ATOM_ATTRS[attr])
    return col

