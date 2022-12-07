import numpy as np
from numpy.typing import ArrayLike

from ptools.namedarray import NamedArrayContainer
from ptools.particlecollection import ParticleCollection

from .testing import assert_array_almost_equal


def test_empty_initialization():
    """Test that the default initialization works."""
    pc = ParticleCollection()
    assert pc.atom_properties is not None
    assert pc.size() == 0


def test_from_attributes():
    """Test initialization from a ``NamedArrayContainer``."""
    attrs = NamedArrayContainer()
    attrs.add_array("x", "xs", [1, 2, 3])
    attrs.add_array("y", "ys", [4, 5, 6])
    attrs.add_array("z", "zs", [7, 8, 9])

    pc = ParticleCollection.from_attributes(attrs)
    assert pc.atom_properties is not None
    assert pc.size() == 3

    assert pc.atom_properties.get("xs") == [1, 2, 3]
    assert pc.atom_properties.get("ys") == [4, 5, 6]
    assert pc.atom_properties.get("zs") == [7, 8, 9]


def test_from_objects():
    """Test initialization from a list of objects."""

    class Dummy:
        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z

    objects = [Dummy(1, 4, 7), Dummy(2, 5, 8), Dummy(3, 6, 9)]

    pc = ParticleCollection.from_objects(objects)
    assert pc.atom_properties is not None
    assert pc.size() == 3

    assert pc.atom_properties.get("xs") == [1, 2, 3]
    assert pc.atom_properties.get("ys") == [4, 5, 6]
    assert pc.atom_properties.get("zs") == [7, 8, 9]


import random


class DummyAtom:
    """Dummy atom class for testing."""


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


def random_amino_acid_name():
    """Returns a random amino acid name."""
    return random.choice(AMINO_ACID_NAMES)


def random_atom_name():
    """Returns a random atom name."""
    return random.choice(["H", "C", "N", "O", "S"])


def random_charge(lower: float = -1, upper: float = 1):
    """Returns a random charge."""
    return random.uniform(lower, upper)


def generate_atoms(size: int = 10) -> list:
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

    class DummyAtom:
        """Dummy atom class for testing."""

        def __init__(self, index, coordinates: ArrayLike):
            self.chain = "A"
            self.index = index
            self.coordinates = np.array(coordinates, dtype=float)
            self.name = random_atom_name()
            self.residue_name = random_amino_acid_name()
            self.charge = random_charge()

    return [DummyAtom(i, [i, i, i]) for i in range(size)]


class TestParticleCollection:
    """Test the ``ParticleCollection`` class."""
    def setup_method(self):
        self.atoms = generate_atoms()

    @property
    def chains(self):
        return [atom.chain for atom in self.atoms]

    @property
    def indices(self):
        return [atom.index for atom in self.atoms]

    @property
    def coordinates(self):
        return np.array([atom.coordinates for atom in self.atoms])

    @property
    def names(self):
        return [atom.name for atom in self.atoms]

    @property
    def residue_names(self):
        return [atom.residue_name for atom in self.atoms]

    @property
    def charges(self):
        return [atom.charge for atom in self.atoms]


    def test_from_atoms(self):
        """Test initialization from a list of AtomAttrs."""
        from ptools.spelling import plural

        pc = ParticleCollection.from_objects(self.atoms)
        assert pc.atom_properties is not None
        assert pc.size() == len(self.atoms)

        for attr in vars(self.atoms[0]):
            name = plural(attr)
            print(name)
            assert pc.atom_properties.get(name) == getattr(self, name)


# def test_get_slice():
#     atoms = ParticleCollection.from_objects(generate_atoms())
#     assert atoms.size() == 3
#     subset = atoms[1:3]
#     assert subset.size() == 2
#     assert subset[0] == atoms[1]
#     assert subset[1] == atoms[2]


# def test_iter():
#     atoms = generate_atoms()
#     col = ParticleCollection.from_objects(atoms)
#     attribute_names = vars(atoms[0])
#     for expected, actual in zip(atoms, col):
#         for name in attribute_names:
#             if name == "coordinates":
#                 assert_array_almost_equal(
#                     getattr(expected, name), getattr(actual, name)
#                 )
#             else:
#                 assert getattr(expected, name) == getattr(actual, name)


# def test_len():
#     """Test that the ``len`` function works."""
#     pc = ParticleCollection()
#     assert len(pc) == 0

#     pc = ParticleCollection.from_objects(generate_atoms())
#     assert len(pc) == 3
