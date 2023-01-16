import random

import numpy as np
from numpy.typing import ArrayLike

import pytest

from ptools.namedarray import NamedArrayContainer
from ptools.particlecollection import Particle, ParticleCollection

from .testing import assert_array_almost_equal


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


class DummyAtom:
    """Dummy atom class for testing."""

    def __init__(self, index, coordinates: ArrayLike):
        self.chain = "A"
        self.index = index
        self.residue_index = index
        self.coordinates = np.array(coordinates, dtype=float)
        self.name = random_atom_name()
        self.residue_name = random_amino_acid_name()
        self.charge = random_charge()


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
    return [DummyAtom(i, [i, i, i]) for i in range(size)]


# ==============================================================================
class TestParticle:
    """Test the ``Particle`` class."""

    def test_equality(self):
        """Tests that two particles are equal if they have the same attributes."""
        col = ParticleCollection.from_objects(generate_atoms(2))
        a1 = col.particles[0]
        a2 = col.particles[0]
        assert a1 == a2

    def test_inequality(self):
        """Tests that two particles are not equal if they have different attributes."""
        col = ParticleCollection.from_objects(generate_atoms(2))
        a1 = col.particles[0]
        a2 = col.particles[1]
        assert a1 != a2

    def test_equality_non_particle(self):
        """Tests that a particle is equal to a non-particle object with the same attributes."""
        class Dummy:
            pass

        col = ParticleCollection.from_objects(generate_atoms(2))
        a1 = col.particles[0]
        a2 = Dummy()
        for attr in col.atom_properties.singular_names():
            setattr(a2, attr, getattr(a1, attr))
        assert a1 == a2

    def test_inequality_non_particle(self):
        """Tests that a particle is not equal to a non-particle object with different attributes."""
        class Dummy:
            pass

        col = ParticleCollection.from_objects(generate_atoms(2))
        a1 = col.particles[0]
        a2 = Dummy()
        for attr in col.atom_properties.singular_names():
            setattr(a2, attr, "X")
        assert a1 != a2

    def test_inequality_missing_attribute(self):
        """Tests that a particle particle equality returns False and does not fail when
        comparing two particles with different attributes (attribute not present)."""
        class Dummy:
            pass

        col = ParticleCollection.from_objects(generate_atoms(2))
        a1 = col.particles[0]
        a2 = Dummy()
        assert a1 != a2

# ==============================================================================

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

    pc = ParticleCollection(attrs)
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
    def residue_indices(self):
        return [atom.residue_index for atom in self.atoms]

    @property
    def charges(self):
        return [atom.charge for atom in self.atoms]

    def test_initialization_from_list_fails(self):
        with pytest.raises(TypeError):
            ParticleCollection(self.atoms)

    def test_from_atoms(self):
        """Test initialization from a list of AtomAttrs.

        Checks that the atoms attributes are correctly stored in the
        ``ParticleCollection`` in a ``NamedArrayContainer`` instance.
        """
        from ptools.spelling import plural

        pc = ParticleCollection.from_objects(self.atoms)
        assert pc.atom_properties is not None
        assert pc.size() == len(self.atoms)

        for attr in vars(self.atoms[0]):
            name = plural(attr)
            assert pc.atom_properties.get(name) == getattr(self, name)

    def test_get_slice(self):
        atoms = ParticleCollection.from_objects(self.atoms)
        assert atoms.size() >= 10
        subset = atoms[1:3]
        assert subset.size() == 2
        assert subset[0] == atoms[1]
        assert subset[1] == atoms[2]


    def test_get_property_array(self):
        """Tests access to the properties of the atoms."""
        atoms = ParticleCollection.from_objects(self.atoms)
        assert atoms.atom_properties.get("chains") == self.chains
        assert atoms.atom_properties.get("indices") == self.indices
        assert atoms.atom_properties.get("coordinates") == self.coordinates
        assert atoms.atom_properties.get("names") == self.names
        assert atoms.atom_properties.get("residue_names") == self.residue_names
        assert atoms.atom_properties.get("residue_indices") == self.residue_indices
        assert atoms.atom_properties.get("charges") == self.charges


    def test_iter(self):
        col = ParticleCollection.from_objects(self.atoms)
        attribute_names = vars(self.atoms[0])
        for expected, actual in zip(self.atoms, col):  # calls __iter__
            for name in attribute_names:
                if name == "coordinates":
                    assert_array_almost_equal(
                        getattr(expected, name), getattr(actual, name)
                    )
                else:
                    assert getattr(expected, name) == getattr(actual, name)

    def test_particle_collection_coordinates(self):
        """Test that the ``coordinates`` property works."""
        pc = ParticleCollection.from_objects(self.atoms)
        assert_array_almost_equal(pc.coordinates, self.coordinates)

    def test_len(self):
        """Test that the ``len`` function works."""
        pc = ParticleCollection()
        assert len(pc) == 0

        pc = ParticleCollection.from_objects(self.atoms)
        assert len(pc) == len(self.atoms)

    def test_contains(self):
        """Test that the ``in`` operator works."""
        pc = ParticleCollection.from_objects(self.atoms)
        for atom in self.atoms:
            assert atom in pc

    def test_getitem_with_int_returns_particle(self):
        """Test that the ``getitem`` operator works with an integer and returns a ``Particle``."""
        pc = ParticleCollection.from_objects(self.atoms)
        for i, atom in enumerate(self.atoms):
            assert isinstance(pc[i], Particle)
            assert pc[i] == atom

    def test_iter_returns_reference(self):
        """Test that the ``in`` operator returns a reference to the original object."""
        pc = ParticleCollection.from_objects(self.atoms)
        for p in pc:
            p.name = "banana"
        expected = ["banana"] * len(self.atoms)
        actual = pc.atom_properties.get("names")
        assert actual == expected

    def test_getitem_with_int_returns_reference(self):
        """Test that the ``getitem`` operator returns a reference to the original object."""
        pc = ParticleCollection.from_objects(self.atoms)
        pc[0].name = "banana"
        assert pc[0].name == "banana"

    def test_copy(self):
        """Test that the ``copy`` method returns a deep copy."""
        pc = ParticleCollection.from_objects(self.atoms)
        pc_copy = pc.copy()
        assert pc_copy is not pc
        assert pc_copy == pc

        pc_copy[0].name = "XXX"
        assert pc_copy[0].name != pc[0].name


def test_add():
    """Test that the ``add`` method works."""
    pc1 = ParticleCollection.from_objects(generate_atoms())
    pc2 = ParticleCollection.from_objects(generate_atoms())

    pc3 = pc1 + pc2

    assert pc3.size() == pc1.size() + pc2.size()
    assert pc3[0] == pc1[0]
    assert pc3[-1] == pc2[-1]
