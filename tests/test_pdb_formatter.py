from ptools.atomattrs import AtomAttrs
from ptools.io import to_pdb

from .generators import generate_atom_attrs

DEFAULT_ATOM_ATTRS = {
    "name": "CA",
    "index": 42,
    "residue_name": "ALA",
    "residue_index": 17,
    "chain": "A",
    "coordinates": (1, 2, 3),
}


def generate_atom(**kwargs) -> AtomAttrs:
    """Creates an ``AtomAttrs`` instance."""
    attrs = DEFAULT_ATOM_ATTRS.copy()
    attrs.update(kwargs)
    return generate_atom_attrs(**attrs)


def test_topdb():
    atom = generate_atom_attrs()
    reference_string = (
        "ATOM     42  CA  ALA A  17       1.000   2.000   3.000  1.00  0.00           X"
    )
    assert to_pdb(atom) == reference_string

    atom.guess_element()
    reference_string = reference_string[:-1] + "C"
    assert to_pdb(atom) == reference_string


def test_topdb_long_atomid():
    atom = generate_atom(index=110000)
    reference_string = (
        "ATOM  1adb0  CA  ALA A  17       1.000   2.000   3.000  1.00  0.00           X"
    )
    assert to_pdb(atom) == reference_string


def test_topdb_long_resid():
    atom = generate_atom(residue_index=11000)
    reference_string = (
        "ATOM     42  CA  ALA A2af8       1.000   2.000   3.000  1.00  0.00           X"
    )
    assert to_pdb(atom) == reference_string


def test_topdb_long_atom_name():
    atom = generate_atom_attrs(name="CA1")
    reference_string = (
        "ATOM     42  CA1 ALA A  17       1.000   2.000   3.000  1.00  0.00           X"
    )
    assert to_pdb(atom) == reference_string


def test_to_pdb_list_of_atoms():
    atoms = [generate_atom() for _ in range(10)]
    result = to_pdb(atoms)
    for line, atom in zip(result.splitlines(), atoms):
        assert line == to_pdb(atom)
