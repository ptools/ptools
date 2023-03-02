"""Tests for the ``ParticleCollection`` class."""

import numpy as np

import pytest
from pytest import approx

import ptools
from ptools.namedarray import NamedArray, NamedArrayContainer
from ptools.particlecollection import Particle, ParticleCollection

from .generators import generate_atoms
from .testing import assert_array_almost_equal, assert_array_equal
from . import TEST_LIGAND


class RandomParticleContainer:
    """Stores a list of dummy particles and provides accessors to their properties."""

    def __init__(self, size: int = 5):
        self.particles = generate_atoms(size)

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

    def number_of_particles(self) -> int:
        return len(self.particles)

    def number_of_residues(self) -> int:
        return len(set(self.residue_indices))

    def number_of_chains(self) -> int:
        return len(set(self.chains))


# == Initialization  ===================================================================


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
    assert pc.size() == expected.number_of_particles()

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


# == __getitem__  ======================================================================


def test_getitem_with_int_returns_a_particle():
    """Test that the ``getitem`` operator works with an integer and returns a ``Particle``."""
    expected = RandomParticleContainer()
    pc = ParticleCollection(expected.particles)
    assert pc.size() == expected.number_of_particles()
    for i, atom in enumerate(expected.particles):
        assert isinstance(pc[i], Particle)
        assert pc[i] == atom


def test_getitem_with_int_returns_reference_to_the_original_object():
    """Test that the ``getitem`` operator returns a reference to the original object."""
    pc = ParticleCollection(generate_atoms())
    pc[0].name = "banana"
    assert pc[0].name == "banana"

    pc[0].coordinates = [12, 12, 12]
    assert pc[0].coordinates == approx([12, 12, 12])


def test_getitem_with_slice():
    expected = RandomParticleContainer(10)
    atoms = ParticleCollection(expected.particles)
    assert expected.number_of_particles() == 10
    assert atoms.size() == expected.number_of_particles()

    # Creates a subset of the atoms.
    subset = atoms[1:3]

    # Checks that the subset has been successfully created.
    assert isinstance(subset, ParticleCollection)
    assert subset.size() == 2
    for i, atom in enumerate(subset):
        assert isinstance(atom, Particle)
        assert atom == expected.particles[i + 1]


def test_subset_with_int_returns_a_reference_to_the_original_object():
    """Test that the ``subset`` method returns a reference to the original object."""
    pc = ParticleCollection(generate_atoms())
    atom = pc[0]
    atom.name = "banana"
    assert pc[0].name == "banana"


def test_subset_with_slice_returns_a_reference_to_the_original_object():
    """Test that the ``subset`` method returns a reference to the original object."""
    pc = ParticleCollection(generate_atoms())
    subset = pc[0:2]
    for atom in subset:
        atom.name = "banana"
    assert pc.names[0:2].tolist() == ["banana", "banana"]


# == __setitem__  ======================================================================


def test_set_property():
    expected = RandomParticleContainer()
    atoms = ParticleCollection(expected.particles)

    expected_names = list(expected.names)
    assert_array_equal(atoms.names, expected_names)

    for i, atom in enumerate(atoms):
        atom.name = atom.name + "random suffix"
        expected_names[i] = expected_names[i] + "random suffix"

    assert_array_equal(atoms.names, expected_names)


def test_set_property_with_inplace_add():
    expected = RandomParticleContainer()
    atoms = ParticleCollection(expected.particles)

    expected_names = list(expected.names)
    assert_array_equal(atoms.names, expected_names)

    for i, atom in enumerate(atoms):
        atom.name += "random suffix"
        expected_names[i] += "random suffix"

    assert_array_equal(atoms.names, expected_names)


def test_set_property_with_getitem():
    expected = RandomParticleContainer()
    atoms = ParticleCollection(expected.particles)

    expected_names = list(expected.names)
    assert_array_equal(atoms.names, expected_names)

    for i in range(len(atoms)):
        atoms[i].name = atoms[i].name + "random suffix"
        expected_names[i] = expected_names[i] + "random suffix"

    assert_array_equal(atoms.names, expected_names)


def test_set_property_with_getitem_and_inplace_add():
    expected = RandomParticleContainer()
    atoms = ParticleCollection(expected.particles)

    expected_names = list(expected.names)
    assert_array_equal(atoms.names, expected_names)

    for i in range(len(atoms)):
        atoms[i].name += "random suffix"
        expected_names[i] += "random suffix"

    assert_array_equal(atoms.names, expected_names)


def test_set_property_using_slice():
    expected = RandomParticleContainer(10)
    atoms = ParticleCollection(expected.particles)
    assert expected.number_of_particles() == 10
    assert atoms.size() == expected.number_of_particles()

    expected_names = list(expected.names)
    assert_array_equal(atoms.names, expected_names)

    atoms.names[1:3] = [name + "random suffix" for name in atoms[1:3].names]
    expected_names[1:3] = [name + "random suffix" for name in expected_names[1:3]]

    assert_array_equal(atoms.names, expected_names)


# == Container methods  ================================================================


def test_size_and_len():
    """Test that the ``len`` and ``size`` functions work."""
    expected = RandomParticleContainer()
    atoms = ParticleCollection(expected.particles)

    assert len(atoms) == expected.number_of_particles()
    assert atoms.size() == expected.number_of_particles()


def test_contains():
    """Test that the ``in`` operator works."""
    expected = RandomParticleContainer()
    atoms = ParticleCollection(expected.particles)
    for atom in expected.particles:
        assert atom in atoms


# == ParticleCollection merging  ========================================================


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


# ======================================================================================


def test_copy():
    """Test that the ``copy`` method returns a deep copy."""
    pc = ParticleCollection(generate_atoms())
    pc_copy = pc.copy()
    assert pc_copy is not pc
    assert pc_copy == pc

    pc_copy[0].name = "XXX"
    assert pc_copy[0].name != pc[0].name


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


# == Selections  ========================================================================


def test_select_atom_type():
    atoms = ptools.read_pdb(TEST_LIGAND)
    sel = atoms.select_atom_type("CA")
    assert len(sel) == 426
    assert sel.names.tolist() == ["CA"] * 426


def test_select_atom_types():
    atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)
    sel = atoms.select_atom_types(["CA", "CB"])
    assert len(sel) == 426 + 64
    assert sorted(sel.names.tolist()) == ["CA"] * 426 + ["CB"] * 64


def test_select_residue_range():
    atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)
    sel = atoms.select_residue_range(10, 20)
    assert len(sel) == 23
    assert all(10 <= atom.residue_index <= 20 for atom in sel)


def test_select_chain():
    atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)
    sel = atoms.select_chain("B")
    assert len(sel) == 974
    assert sel.chains.tolist() == ["B"] * 974


# == Grouping  ==========================================================================


def test_groupby():
    """Test that the ``groupby`` method works.

    Very basic test, just makes sure that grouping by the residue name attribute
    returns the 20 groups (very dependant on the input PDB file).
    """
    atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)
    grouped = atoms.groupby(key=lambda atom: atom.residue_name)
    assert len(grouped) == 20


# == Iteration methods  =================================================================


def test_iter():
    expected = RandomParticleContainer()
    atoms = ParticleCollection(expected.particles)
    assert atoms.size() == expected.number_of_particles()

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


def test_iter_residues():
    """Test that the ``iter_residues`` method works."""
    expected = RandomParticleContainer()
    pc = ParticleCollection(expected.particles)

    iterator = pc.iter_residues()
    # assert


# == Add/remove atom properties  ========================================================


def test_add_atom_property():
    """Test that the ``add_property`` method works."""
    pc = ParticleCollection(generate_atoms())
    bananas = [0] * pc.size()
    pc.add_atom_property("banana", "bananas", bananas)
    assert hasattr(pc, "bananas")
    assert pc.atom_properties.get("bananas") == bananas


def test_add_atom_property_fails_if_property_already_exists():
    """Test that the ``add_property`` method fails if the property already exists."""
    pc = ParticleCollection(generate_atoms())
    bananas = [0] * pc.size()
    pc.add_atom_property("banana", "bananas", bananas)
    with pytest.raises(KeyError):
        pc.add_atom_property("banana", "bananas", bananas)


def test_add_atom_property_fails_when_propery_has_wrong_dimensions():
    """Test that the ``add_property`` method fails if the property has the wrong dimensions."""
    pc = ParticleCollection(generate_atoms())
    bananas = [0] * (pc.size() - 1)
    with pytest.raises(ValueError):
        pc.add_atom_property("banana", "bananas", bananas)


def test_remove_atom_property():
    """Test that the ``remove_property`` method works."""
    pc = ParticleCollection(generate_atoms())
    assert hasattr(pc, "names")
    pc.remove_atom_property("names")
    assert not hasattr(pc, "names")
