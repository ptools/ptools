
"""pyptools.atom - Defines classes and function that handle atom and
atom groups."""

import copy

import numpy

from pyptools.spatial import SpatialObject, coord3d


class BaseAtom(SpatialObject):
    """Base class for an Atom.

    Args:
        index (int):
        name (str):
        resname (str):
        chain (str):
        resid (int):
        charge (float):
        coords (numpy.ndarray):
        orig (BaseAtom): initialize from other BaseAtom (copy constructor)
    """

    def __init__(self, index=0, name='XXX', resname='XXX', chain='X',
                 resid=0, charge=0.0, coords=(0, 0, 0), orig=None):
        if orig is not None:
            self._init_copy(orig)
        else:
            self._init_non_copy(index, name, resname, chain, resid, charge,
                                coords)

    def _init_non_copy(self, index, name, resname, chain, resid, charge,
                       coords):
        super().__init__(coords)
        self.name = name
        self.resname = resname
        self.chain = chain
        self.index = index
        self.resid = resid
        self.charge = charge

    def _init_copy(self, other):
        super().__init__(other.coords)
        self.name = other.name
        self.resname = other.resname
        self.chain = other.chain
        self.index = other.index
        self.resid = other.resid
        self.charge = other.charge

    def copy(self):
        """Return a copy of the current atom."""
        return copy.copy(self)


class Atom(BaseAtom):
    """Atom that belongs to a group of atoms.

    Its coordinates are a weak reference to the AtomCollection coordinate array.
    Atom knows its positions in the AtomCollection and can therefore find
    its coordinates into the AtomCollection coordinate array.

    An Atom initialization can only be done from a BaseAtom.

    Args:
        serial (int): atom index in the collection (can be different from
            `Atom.index` attribute).
        atom (BaseAtom): atom properties stored in BaseAtom.
        collection (AtomCollection): reference to the collection the atom
            belongs to.

    """
    def __init__(self, atom, serial=0, collection=None):
        super().__init__(orig=atom)
        self.serial = serial
        self.collection = collection



class AtomCollection(object):
    """Group of atoms.

    For better performances, atom coordinates are stored into a numpy array.

    Args:
        atoms (list[BaseAtom]): list of atoms
    """
    def __init__(self, atoms):
        self.atoms = [Atom(atom, serial, self)
                      for serial, atom in enumerate(atoms)]

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)






