
"""pyptools.atom - Defines classes and function that handle atom and
atom groups."""

import copy
import math

import numpy as np

from pyptools.spatial import SpatialObject, coord3d


# The Protein Data Bank format for atom coordinates
PDB_FMT = "%-6s%5s %4s%c%-4s%c%4s%s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s"


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

    def topdb(self):
        def first_alpha(s):
            for c in s:
                if c.isalpha():
                    return c
            return ""

        rec = "ATOM"
        indexbuf = "*****"
        residbuf = "****"
        altlocchar = " "
        insertion = " "
        chain = " " if not self.chain else self.chain
        occ = 1.0
        b = 0.0
        element = first_alpha(self.name)

        if self.index < 100000:
            indexbuf = "{0:5d}".format(self.index)
        elif self.index < 1048576:
            indexbuf = "{0:05x}".format(self.index)

        if self.resid < 10000:
            residbuf = "{0:4d}".format(self.resid)
        elif self.resid < 65536:
            residbuf = "{0:04x}".format(self.resid)

        namebuf = self.name.center(4)
        if len(self.name) > 2:
            namebuf = "{0:>4s}".format(self.name)

        x, y, z = self.coords

        return PDB_FMT % (rec, indexbuf, namebuf, altlocchar,
                          self.resname, chain[0], residbuf, insertion[0],
                          x, y, z, occ, b, "", element)


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
        if self.atoms:
            coords = [atom._coords for atom in self.atoms]
        else:
            coords = np.zeros((0, 3))
        super().__init__(coords)

    def copy(self):
        """Return a copy of the current collection."""
        return self.__class__(self.atoms)

    def __len__(self):
        """Get the number of atoms in the collection."""
        return len(self.atoms)

    def size(self):
        """Get the number of atoms in the collection.

        Alias for len(AtomCollection).
        """
        return len(self)

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
        return self.centroid()

    def get_radius_of_gyration(self):
        """Return the isometric radius of gyration (atom mass is not taken
        into account)."""
        centered = self.coords - self.get_center()
        rgyr2 = np.sum(centered ** 2) / len(self)
        return math.sqrt(rgyr2)

    def topdb(self):
        """Return a string representing the AtomCollection in PDB format."""
        return "\n".join(atom.topdb() for atom in self)

    def writepdb(self, path):
        """Write the AtomCollection to a PDB formatted file."""
        with open(path, "wt") as f:
            print(self.topdb(), file=f)
