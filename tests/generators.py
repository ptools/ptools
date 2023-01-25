"""Generators for testing."""

from dataclasses import dataclass, field
import random
from typing import Optional

import numpy as np
from numpy.typing import ArrayLike

from ptools.atomattrs import AtomAttrs
from ptools.particlecollection import ParticleCollection


AMINO_ACID_NAMES = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]


def generate_atom_attrs(**kwargs) -> AtomAttrs:
    """Creates an ``AtomAttrs`` instance."""
    attrs = {
        "name": "CA",
        "index": 42,
        "residue_name": "ALA",
        "residue_index": 17,
        "chain": "A",
        "charge": 2.0,
        "coordinates": (1, 2, 3),
    }
    attrs.update(kwargs)
    atom = AtomAttrs()
    for attr, value in attrs.items():
        setattr(atom, attr, value)
    return atom


class Atom:
    """Dummy atom class for testing."""

    def __init__(self, index, coordinates: ArrayLike):
        self.chain = "A"
        self.index = index
        self.residue_index = index
        self.coordinates = np.array(coordinates, dtype=float)
        self.name = random_atom_name()
        self.residue_name = random_amino_acid_name()
        self.charge = random_charge()


def generate_atoms(size: int = 10, names: Optional[list[str]] = None) -> list[Atom]:
    """Creates a dummy atom collection composed of `size` atoms.

    The atoms have those properties:

    - coordinates are [(0, 0, 0), (1, 1, 1), ..., (size - 1, size - 1, size - 1)],
    - atom names are randomly chosen from the list of atom names,
    - atom indices are [1, 2, ..., size],
    - residue names are randomly chosen from the list of amino acid names,
    - residue indices are [1, 2, ..., size],
    - chain is "A",
    - charge is a random float between -1 and 1,
    """
    atoms = [Atom(i, [i, i, i]) for i in range(size)]
    if names is not None:
        for atom, name in zip(atoms, names):
            atom.name = name
    return atoms


def generate_particlecollection(**kwargs) -> ParticleCollection:
    """Creates a dummy particle collection."""
    return ParticleCollection(generate_atoms(**kwargs))


def random_amino_acid_name():
    """Returns a random amino acid name."""
    return random.choice(AMINO_ACID_NAMES)


def random_atom_name():
    """Returns a random atom name."""
    return random.choice(["H", "C", "N", "O", "S"])


def random_charge(lower: float = -1, upper: float = 1):
    """Returns a random charge."""
    return random.uniform(lower, upper)


@dataclass
class Balloon:
    """Dummy object with coordinates."""

    coordinates: np.ndarray = field(default_factory=lambda: np.zeros((5, 3)))


def generate_balloon(coordinates: Optional[np.typing.ArrayLike] = None) -> Balloon:
    """Returns a ``Balloon`` with given coordinates"""
    if coordinates is None:
        return Balloon()
    return Balloon(np.asarray(coordinates))
