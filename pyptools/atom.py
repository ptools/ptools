
"""pyptools.atom - Defines classes and function that handle atom and
atom groups."""

import copy
import math

import numpy

from pyptools.spatial import SpatialObject, coord3d, translate


class BaseAtom(SpatialObject):
    """Base class for an Atom.

    Args:
        index (int): index (read from topology file)
        name (str): atom name
        resname (str): residue name
        chain (str): chain identifier
        resid (int): residue number (read from topology file)
        charge (float): charge
        coords (numpy.ndarray): cartesian coordinates
        meta (dict[str, ()]): metadata dictionnary
        orig (BaseAtom): initialize from other BaseAtom (copy constructor)
    """

    def __init__(self, index=0, name='XXX', resname='XXX', chain='X',
                 resid=0, charge=0.0, coords=(0, 0, 0), meta={}, orig=None):
        if orig is not None:
            self._init_copy(orig)
        else:
            self._init_non_copy(index, name, resname, chain, resid, charge,
                                coords, meta)

    def _init_non_copy(self, index, name, resname, chain, resid, charge,
                       coords, meta):
        super().__init__(coords)
        self.name = name
        self.resname = resname
        self.chain = chain
        self.index = index
        self.resid = resid
        self.charge = charge
        self.meta = meta

    def _init_copy(self, other):
        super().__init__(other.coords)
        self.name = other.name
        self.resname = other.resname
        self.chain = other.chain
        self.index = other.index
        self.resid = other.resid
        self.charge = other.charge
        self.meta = other.meta

    def copy(self):
        """Return a copy of the current atom."""
        return copy.copy(self)


class Atom(BaseAtom):
    """Atom that belongs to a group of atoms.

    Its coordinates are a weak reference to the AtomCollection coordinate
    array. Atom knows its positions in the AtomCollection and can therefore
    find its coordinates into the AtomCollection coordinate array.

    An Atom initialization can only be done from a BaseAtom.

    Args:
        serial (int): atom index in the collection (can be different from
            `Atom.index` attribute).
        atom (BaseAtom): atom properties stored in BaseAtom.
        collection (AtomCollection): reference to the collection the atom
            belongs to.

    """
    def __init__(self, atom, serial, collection):
        super().__init__(orig=atom)
        self.serial = serial
        self.collection = collection

    @property
    def coords(self):
        """Get atom cartesian coordinates."""
        return self.collection.coords[self.serial].copy()

    @coords.setter
    def coords(self, pos):
        self.collection.coords[self.serial] = coord3d(pos)


class AtomCollection(SpatialObject):
    """Group of atoms.

    For better performances, atom coordinates are stored into a numpy array.

    Args:
        atoms (list[BaseAtom]): list of atoms
    """
    def __init__(self, atoms):
        self.atoms = [Atom(atom, serial, self)
                      for serial, atom in enumerate(atoms)]
        super().__init__([atom._coords for atom in self.atoms])

    def copy(self):
        """Return a copy of the current collection."""
        return self.__class__(self.atoms)

    def __len__(self):
        """Get the number of atoms in the collection."""
        return len(self.atoms)

    def __iter__(self):
        """Iterate over the collection atoms."""
        return iter(self.atoms)

    def __getitem__(self, serial):
        """Access an atom by its serial number (which the internal index
        starting at 0)."""
        return self.atoms[serial]

    def get_center(self):
        """Return the isobarycenter (geometric center) of a collection of
        atoms."""
        return sum(self.coords) / len(self)

    def get_radius_of_gyration(self):
        """Return the isometric radius of gyration (atom mass is not taken
        into account)."""
        centered = self.coords - self.get_center()
        rgyr2 = numpy.sum(centered ** 2) / len(self)
        return math.sqrt(rgyr2)
