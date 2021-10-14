"""Creates dummy ptools structures."""

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


def dummy_atom():
    """Creates a dummy atom."""
    atom = BaseAtom()
    for attr, value in DUMMY_ATOM_ATTRS.items():
        setattr(atom, attr, value)
    return atom


def dummy_atomcollection(n_atoms:int = 10):
    """Creates a dummy atom collection composed of `n_atoms` atoms.

    Atom coordinates are [(0, 0, 0), (1, 1, 1), ..., (9, 9, 9)].
    """
    return AtomCollection([BaseAtom(coords=(i, i, i)) for i in range(n_atoms)])
