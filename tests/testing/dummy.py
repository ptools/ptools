"""Creates dummy ptools structures."""

from ptools.atom import Atom, BaseAtom, AtomCollection
from ptools.mcop import Mcop, McopRigid

DUMMY_ATOM_ATTRS = {
    "name": "CA",
    "resname": "ALA",
    "chain": "A",
    "index": 42,
    "resid": 17,
    "charge": 2.0,
    "coords": (1, 2, 3),
}


def dummy_atom() -> BaseAtom:
    """Creates a dummy atom."""
    atom = BaseAtom()
    for attr, value in DUMMY_ATOM_ATTRS.items():
        setattr(atom, attr, value)
    return atom


def dummy_atomcollection(n_atoms:int = 10) -> AtomCollection:
    """Creates a dummy atom collection composed of `n_atoms` atoms.

    Atom coordinates are [(0, 0, 0), (1, 1, 1), ..., (9, 9, 9)].
    """
    return AtomCollection([BaseAtom(coords=(i, i, i)) for i in range(n_atoms)])


def dummy_mcop(n_copies: int = 10) -> Mcop:
    """Creates a dummpy Mcop composed of `n_copies` dummy AtomCollection instances."""
    mcop = Mcop()
    for _ in range(n_copies):
        mcop.add_copy(dummy_atomcollection())
    return mcop


def dummy_mcop_rigid(n_regions: int = 10) -> McopRigid:
    """Creates a dummpy McopRigid composed of `n_regions` that are dummy Mcop instances."""
    rigid = McopRigid()
    rigid.core = dummy_atomcollection()
    for _ in range(n_regions):
        rigid.add_region(dummy_mcop())
    return rigid