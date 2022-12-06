
from ptools.atomattrs import AtomAttrs
from ptools.namedarray import NamedArrayContainer
from ptools.particlecollection import ParticleCollection


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



def generate_atoms():
    from ptools.atomattrs import AtomAttrs
    atoms = [
        AtomAttrs(index=1, name="H", residue_name="ALA", coordinates=(1, 2, 3)),
        AtomAttrs(index=2, name="C", residue_name="VAL", coordinates=(4, 5, 6)),
        AtomAttrs(index=3, name="N", residue_name="HIS", coordinates=(7, 8, 9)),
    ]
    return atoms


def test_from_atoms():
    """Test initialization from a list of AtomAttrs."""
    atoms = generate_atoms()
    pc = ParticleCollection.from_objects(atoms)
    assert pc.atom_properties is not None
    assert pc.size() == 3

    expected_coordinates = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ]
    assert pc.atom_properties.get("coordinates") == expected_coordinates

    assert pc.atom_properties.get("indices") == [1, 2, 3]
    assert pc.atom_properties.get("names") == ["H", "C", "N"]
    assert pc.atom_properties.get("residue_names") == ["ALA", "VAL", "HIS"]


def test_get_slice():
    atoms = ParticleCollection.from_objects(generate_atoms())
    assert atoms.size() == 3
    subset = atoms[1:3]
    assert subset.size() == 2

    for a in atoms:
        print(a)

    assert 1 == 2
