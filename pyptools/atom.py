
"""pyptools.atom - Defines classes and function that handle atom and
atom groups."""

import copy
import math

import numpy as np

from pyptools.spatial import SpatialObject, coord3d
from . import tables


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
                 resid=0, charge=0.0, coords=(0, 0, 0), meta={},
                 orig=None):
        if orig is not None:
            self._init_copy(orig)
        else:
            self._init_non_copy(index, name, resname, chain, resid, charge,
                                coords, meta)

    def _init_non_copy(self, index, name, resname, chain, resid, charge,
                       coords, meta):
        super().__init__(coords)
        self._name = name
        self.resname = resname
        self.chain = chain
        self.index = index
        self.resid = resid
        self.charge = charge
        self.meta = meta
        self.element = guess_atom_element(self.name)
        self.mass = guess_atom_mass(self.element)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, s):
        self._name = s
        self.element = guess_atom_element(self.name)
        self.mass = guess_atom_mass(self.element)

    def _init_copy(self, other):
        super().__init__(other.coords)
        self._name = other._name
        self.resname = other.resname
        self.chain = other.chain
        self.index = other.index
        self.resid = other.resid
        self.charge = other.charge
        self.meta = other.meta
        self.element = other.element
        self.mass = other.mass

    def copy(self):
        """Returns a copy of the current atom."""
        return copy.copy(self)

    def topdb(self):
        rec = "ATOM"
        indexbuf = "*****"
        residbuf = "****"
        altlocchar = " "
        insertion = " "
        chain = " " if not self.chain else self.chain
        occ = 1.0
        b = 0.0
        element = self.element

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
        """Gets atom cartesian coordinates."""
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
    def __init__(self, atoms=[]):
        self.atoms = [Atom(atom, serial, self)
                      for serial, atom in enumerate(atoms)]
        if self.atoms:
            coords = [atom._coords for atom in self.atoms]
        else:
            coords = np.zeros((0, 3))
        super().__init__(coords)

    def masses(self):
        """Returns the array of atom masses."""
        return np.array([atom.mass for atom in self.atoms])

    def copy(self):
        """Returns a copy of the current collection."""
        return self.__class__(self.atoms)

    def __len__(self):
        """Gets the number of atoms in the collection."""
        return len(self.atoms)

    def size(self):
        """Gets the number of atoms in the collection.

        Alias for len(AtomCollection).
        """
        return len(self)

    def __iter__(self):
        """Iterate over the collection atoms."""
        return iter(self.atoms)

    def __getitem__(self, serial):
        """Accesses an atom by its serial number (which the internal index
        starting at 0)."""
        return self.atoms[serial]

    def __add__(self, other):
        """Concatenate two RigidBody instances."""
        output = self.copy()
        for atom in other:
            output.atoms.append(atom.copy())
            output._coords = np.vstack((output._coords, atom._coords))
        return output

    def center(self):
        """Returns the isobarycenter (geometric center) of a collection of
        atoms."""
        return self.centroid()

    def center_of_mass(self):
        """Return the center of mass (barycenter)."""
        weights = self.masses().reshape(-1, 1)
        return (self.coords * weights).sum(axis=0)  / weights.sum()

    def moment_of_inertia(self):
        masses = self.masses()
        com = self.center_of_mass()
        pos = self.coords - com

        tens = np.zeros((3, 3))

        tens[0][0] = (masses * (pos[:, 1] ** 2 + pos[:, 2] ** 2)).sum()
        tens[0][1] = - (masses * pos[:, 0] * pos[:, 1]).sum()
        tens[0][2] = - (masses * pos[:, 0] * pos[:, 2]).sum()

        tens[1][0] = tens[0][1]
        tens[1][1] = (masses * (pos[:, 0] ** 2 + pos[:, 2] ** 2)).sum()
        tens[1][2] = - (masses * pos[:, 1] * pos[:, 2]).sum()

        tens[2][0] = tens[0][2]
        tens[2][1] = tens[1][2]
        tens[2][2] = (masses * (pos[:, 0] ** 2 + pos[:, 1] ** 2)).sum()

        return tens

    def principal_axes(self):
        """Calculate the principal axes."""
        val, vec = np.linalg.eig(self.moment_of_inertia())
        indices = np.argsort(val)[::-1]
        return vec[:, indices].T

    def radius_of_gyration(self):
        """Returns the isometric radius of gyration (atom mass is not taken
        into account)."""
        centered = self.coords - self.center()
        rgyr2 = np.sum(centered ** 2) / len(self)
        return math.sqrt(rgyr2)

    def topdb(self):
        """Returns a string representing the AtomCollection in PDB format."""
        return "\n".join(atom.topdb() for atom in self)

    def writepdb(self, path):
        """Writes the AtomCollection to a PDB formatted file."""
        with open(path, "wt") as f:
            print(self.topdb(), file=f)

    def set_chain(self, chain):
        """Sets all atom chain property."""
        for atom in self.atoms:
            atom.chain = chain


def guess_atom_element(atom_name):
    """Returns the atom element based on its name.

    Basically returns the first non-numeric character in atom_name.
    """
    for c in atom_name:
        if c.isalpha():
            return c
    return "X"


def guess_atom_mass(element):
    """Returns the atom mass based on the element name."""
    return tables.masses.get(element, 1.0)
