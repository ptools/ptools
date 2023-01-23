"""Tests for the ``ParticleCollection`` class."""

import numpy as np

from pytest import approx

from ptools.namedarray import NamedArray, NamedArrayContainer
from ptools.particlecollection import Particle, ParticleCollection

from .generators import generate_atoms
from .testing import assert_array_almost_equal


class RandomParticleContainer:
    """Stores a list of dummy particles and provides accessors to their properties."""

    def __init__(self, size: int = 5):
        self.particles = generate_atoms(size)

    def size(self) -> int:
        return len(self.particles)

    @property
    def chains(self):
        return [atom.chain for atom in self.particles]

    @property
    def indices(self):
        return [atom.index for atom in self.particles]

    @property
    def coordinates(self):
        return np.array([atom.coordinates for atom in self.particles])

    @property
    def names(self):
        return [atom.name for atom in self.particles]

    @property
    def residue_names(self):
        return [atom.residue_name for atom in self.particles]

    @property
    def residue_indices(self):
        return [atom.residue_index for atom in self.particles]

    @property
    def charges(self):
        return [atom.charge for atom in self.particles]

    def properties(self) -> list[str]:
        """Returns the list of particles property names."""
        return [
            "chains",
            "indices",
            "coordinates",
            "names",
            "residue_names",
            "residue_indices",
            "charges",
        ]



def test_initialization_from_empty_list():
    """Test that the default initialization works."""
    pc = ParticleCollection()
    assert pc.atom_properties is not None
    assert pc.size() == 0


def test_initialization_from_list_of_particles():
    """Test initialization from a list of dummy atoms."""
    expected = RandomParticleContainer()
    pc = ParticleCollection(expected.particles)

    assert pc.atom_properties is not None
    assert pc.size() == expected.size()

    for prop in expected.properties():
        assert pc.atom_properties.get(prop) == getattr(expected, prop)


def test_initialization_from_properties():
    """Test initialization from a ``NamedArrayContainer``."""
    attrs = NamedArrayContainer()
    attrs.add_array("x", "xs", [1, 2, 3])
    attrs.add_array("y", "ys", [4, 5, 6])
    attrs.add_array("z", "zs", [7, 8, 9])

    pc = ParticleCollection.from_properties(attrs)
    assert pc.atom_properties is not None
    assert pc.size() == 3

    assert pc.atom_properties.get("xs") == [1, 2, 3]
    assert pc.atom_properties.get("ys") == [4, 5, 6]
    assert pc.atom_properties.get("zs") == [7, 8, 9]


def test_get_slice():
    expected = RandomParticleContainer(10)
    atoms = ParticleCollection(expected.particles)

    assert expected.size() == 10
    assert atoms.size() == expected.size()

    # Creates a subset of the atoms.
    subset = atoms[1:3]

    # Checks that the subset is a new object.
    assert subset.size() == 2
    assert subset[0] == atoms[1]
    assert subset[1] == atoms[2]


# def test_get_slice_returns_a_reference_to_the_original_collection():
#     expected = RandomParticleContainer(10)
#     atoms = ParticleCollection(expected.particles)

#     # Creates a subset of the atoms.
#     subset = atoms[1:3]

#     # Modifies the subset and checks the original object is changed as well.
#     subset_expected_names = [a.name for a in subset]
#     assert subset_expected_names == expected.names[1:3]

#     subset[0].name = "new random name"
#     assert subset[0].name == "new random name"
#     assert atoms[1].name == "new random name"


def test_size_and_len():
    """Test that the ``len`` and ``size`` functions work."""
    expected = RandomParticleContainer()
    atoms = ParticleCollection(expected.particles)

    assert len(atoms) == expected.size()
    assert atoms.size() == expected.size()


def test_contains():
    """Test that the ``in`` operator works."""
    expected = RandomParticleContainer()
    atoms = ParticleCollection(expected.particles)
    for atom in expected.particles:
        assert atom in atoms


def test_getitem_with_int_returns_a_particle():
    """Test that the ``getitem`` operator works with an integer and returns a ``Particle``."""
    expected = RandomParticleContainer()
    pc = ParticleCollection(expected.particles)
    for i, atom in enumerate(expected.particles):
        assert isinstance(pc[i], Particle)
        assert pc[i] == atom


def test_iter():
    expected = RandomParticleContainer()
    atoms = ParticleCollection(expected.particles)
    assert atoms.size() == expected.size()

    for expected_atom, actual_atom in zip(expected.particles, atoms):  # calls __iter__
        assert expected_atom == actual_atom


def test_iter_returs_a_reference_to_original_object():
    """Test that the ``in`` operator returns a reference to the original object."""
    expected = RandomParticleContainer()
    pc = ParticleCollection(expected.particles)

    for expected_atom, actual_atom in zip(expected.particles, pc):
        assert expected_atom == actual_atom

    for actual_atom in pc:
        actual_atom.name = "banana"

    assert pc.atom_properties.get("names") == ["banana"] * len(expected.particles)


def test_iter_returns_a_reference():
    """Test that the ``in`` operator returns a reference to the original object."""
    expected = RandomParticleContainer()
    pc = ParticleCollection(expected.particles)
    for p in pc:
        p.name = "banana"
    expected = ["banana"] * len(expected.particles)
    actual = pc.atom_properties.get("names")
    assert actual == expected


def test_getitem_with_int_returns_reference_to_the_original_object():
    """Test that the ``getitem`` operator returns a reference to the original object."""
    pc = ParticleCollection(generate_atoms())
    pc[0].name = "banana"
    assert pc[0].name == "banana"

    pc[0].coordinates = pc[0].coordinates + 12
    assert pc[0].coordinates == approx([12, 12, 12])


def test_copy():
    """Test that the ``copy`` method returns a deep copy."""
    pc = ParticleCollection(generate_atoms())
    pc_copy = pc.copy()
    assert pc_copy is not pc
    assert pc_copy == pc

    pc_copy[0].name = "XXX"
    assert pc_copy[0].name != pc[0].name


def test_add():
    """Test that the ``+`` operator works."""
    pc1 = ParticleCollection(generate_atoms())
    pc2 = ParticleCollection(generate_atoms())

    pc3 = pc1 + pc2

    assert pc3.size() == pc1.size() + pc2.size()
    assert pc3[0] == pc1[0]
    assert pc3[-1] == pc2[-1]


def test_add_makes_copies():
    """Test that the ``+`` operator makes copies of the particles."""
    pc1 = ParticleCollection(generate_atoms())
    pc2 = ParticleCollection(generate_atoms())

    pc3 = pc1 + pc2

    pc3[0].name = pc1[0].name + "random string"
    assert pc3[0].name != pc1[0].name


def test_inplace_add():
    """Test that the ``+=`` operator works."""
    pc1 = ParticleCollection(generate_atoms(size=10))
    pc2 = ParticleCollection(generate_atoms())

    pc1 += pc2

    assert pc1.size() == pc2.size() * 2
    assert pc1[10] == pc2[0]
    assert pc1[-1] == pc2[-1]


def test_masses():
    """Test that the ``guess_masses`` method."""
    from ptools.tables import atomic_masses
    pc = ParticleCollection(generate_atoms())
    pc.guess_masses()

    expected = [atomic_masses[a.name] for a in pc]
    assert pc.atom_properties.get("masses") == expected


def test_set_property():
    """Test that the ``set_property`` method works."""
    pc = ParticleCollection(generate_atoms())
    assert pc.atom_properties.get("indices") == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    expected = [42] * pc.size()
    pc.atom_properties.set("indices", expected)
    assert isinstance(pc.atom_properties.get("indices"), NamedArray)
    assert pc.atom_properties.get("indices") == expected

