
from ptools.io import to_pdb

from .testing.dummy import generate_dummy_atom, generate_dummy_atomcollection

def test_topdb():
    atom = generate_dummy_atom()
    reference_string = (
        "ATOM     42  CA  ALA A  17       "
        "1.000   2.000   3.000  1.00  0.00           "
        "C"
    )
    assert to_pdb(atom) == reference_string


def test_topdb_long_atomid():
    atom = generate_dummy_atom()
    atom.index = 110000
    reference_string = (
        "ATOM  1adb0  CA  ALA A  17       "
        "1.000   2.000   3.000  1.00  0.00           "
        "C"
    )
    assert to_pdb(atom) == reference_string


def test_topdb_long_resid():
    atom = generate_dummy_atom()
    atom.residue_index = 11000
    reference_string = (
        "ATOM     42  CA  ALA A2af8       "
        "1.000   2.000   3.000  1.00  0.00           "
        "C"
    )
    assert to_pdb(atom) == reference_string


def test_topdb_long_atom_name():
    atom = generate_dummy_atom()
    atom.name = "CA1"
    reference_string = (
        "ATOM     42  CA1 ALA A  17       "
        "1.000   2.000   3.000  1.00  0.00           "
        "C"
    )
    assert to_pdb(atom) == reference_string


def test_to_pdb():
    atoms = generate_dummy_atomcollection()
    result = to_pdb(atoms)
    for line, atom in zip(result.splitlines(), atoms):
        assert line == to_pdb(atom)

