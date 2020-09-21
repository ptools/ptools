
"""ptools.atom - Defines classes and function that handle atom and
atom groups."""

import copy
import math

import numpy as np

import ptools
from ptools.spatial import SpatialObject
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
                 resid=0, charge=0.0, coords=(0, 0, 0), meta=None,
                 orig=None):
        super().__init__(coords)
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

    @property
    def name(self):
        """Gets/Sets atom's name."""
        return self._name

    @name.setter
    def name(self, name):
        self._name = name
        self.element = guess_atom_element(self.name)

    def _init_copy(self, other):
        super().__init__(other.coords)
        self._name = other.name
        self.resname = other.resname
        self.chain = other.chain
        self.index = other.index
        self.resid = other.resid
        self.charge = other.charge
        self.meta = other.meta
        self.element = other.element

    def copy(self):
        """Returns a copy of the current atom."""
        obj = copy.copy(self)
        obj.coords = obj.coords.copy()
        return obj

    def topdb(self):
        """Returns the atom's description in PDB format."""
        rec = "ATOM"
        indexbuf = "*****"
        residbuf = "****"
        altlocchar = " "
        insertion = " "
        chain = " " if not self.chain else self.chain
        occ = 1.0
        bfactor = 0.0
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
                          x, y, z, occ, bfactor, "", element)


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
        self.serial = serial
        self.collection = collection
        super().__init__(orig=atom)

    @property
    def coords(self):
        """Gets atom cartesian coordinates."""
        return self.collection.coords[self.serial].copy()

    @coords.setter
    def coords(self, pos):
        self.collection.coords[self.serial] = np.array(pos)

    @property
    def mass(self):
        return self.collection.masses[self.serial]

    @mass.setter
    def mass(self, mass):
        self.collection.masses[self.serial] = mass

    # Override BaseAtom name setter to manually handle masses
    # (name getter override is mandatory).
    @property
    def name(self):
        """Gets/Sets atom's name."""
        return self._name

    @name.setter
    def name(self, name):
        self._name = name
        self.element = guess_atom_element(self.name)
        self.collection.masses[self.serial] = guess_atom_mass(self.element)


class AtomCollection(SpatialObject):
    """Group of atoms.

    For better performances, atom coordinates are stored into a numpy array.

    Args:
        atoms (list[BaseAtom]): list of atoms
    """
    def __init__(self, atoms=None):
        if atoms is None:
            atoms = []
        self.atoms = [Atom(atom, serial, self)
                        for serial, atom in enumerate(atoms)]
        if self.atoms:
            coords = [atom._coords for atom in self.atoms]
        else:
            coords = np.zeros((0, 3))
        super().__init__(coords)
        self.masses = np.zeros(len(self.atoms))
        self.guess_masses()

    def __len__(self):
        """Gets the number of atoms in the collection."""
        return len(self.atoms)

    def __repr__(self):
        """String representation."""
        modulename = self.__module__
        classname = self.__class__.__name__
        return f"<{modulename}.{classname} with {len(self)} atoms>"

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
        other = other.copy()
        output.atoms += [atom for atom in other]
        output.coords = np.concatenate((output.coords, other.coords), axis=0)
        return output

    def guess_masses(self):
        """Guess atom masses and store them."""
        self.masses = np.array([guess_atom_mass(atom.element) for atom in self])

    def copy(self):
        """Returns a copy of the current collection."""
        return self.__class__(self.atoms)

    def size(self):
        """Gets the number of atoms in the collection.

        Alias for len(AtomCollection).
        """
        return len(self)

    def center(self):
        """Returns the isobarycenter (geometric center) of a collection of
        atoms."""
        return self.centroid()

    def center_of_mass(self):
        """Return the center of mass (barycenter)."""
        return ptools.spatial.center_of_mass(self.coords, self.masses)

    def tensor_of_inertia(self, method="accurate"):
        """Returns the inertia tensors of a set of atoms.

        Args:
            method (str): "fast" or "accurate"
                The "fast" method does not take into account atom masses.
        """
        if method == "fast":
            return ptools.spatial.tensor_of_inertia(self.coords, None, method)
        return ptools.spatial.tensor_of_inertia(self.coords, self.masses, method)

    def principal_axes(self, sort=True, method="accurate"):
        """Returns an AtomCollection principal axes.

        Args:
            sort (bool): sort axes by importance
            method (str): "fast" or "accurate"
                The "fast" method does not take into account atom masses when
                calculating the tensor of inertia, resulting in a less
                accurate result.
                The "fast" method is probably sufficient to calculate axes of
                inertia.
        """
        return ptools.spatial.principal_axes(self.tensor_of_inertia(method), sort)

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

    def select_atom_type(self, atom_type):
        """Returns a sub-collection made of atoms with desired atom type."""
        return self.__class__(atoms=[atom for atom in self if atom.name == atom_type])

    def select_residue_range(self, start, end):
        """Returns a sub-collection made of atoms with desired which residue is within the range."""
        return self.__class__(atoms=[atom for atom in self if start <= atom.resid <= end])



def guess_atom_element(atom_name):
    """Returns the atom element based on its name.

    Basically returns the first non-numeric character in atom_name.
    """
    for char in atom_name:
        if char.isalpha():
            return char
    return "X"


def guess_atom_mass(element):
    """Returns the atom mass based on the element name."""
    return tables.masses.get(element, 1.0)
