"""test_rigidbody - Tests for `ptools.rigidbody module."""

from ptools.rigidbody import RigidBody
from .testing.io import mk_pdb_10_atoms
from .generators import generate_atoms


def generate_rigid_from_pdb():
    with mk_pdb_10_atoms() as temporary_file:
        rigid = RigidBody.from_pdb(temporary_file.name)
    return rigid


def test_initialization_from_list_of_atoms():
    rigid = RigidBody(generate_atoms(10))
    assert len(rigid) == 10


def test_initialization_from_pdb():
    with mk_pdb_10_atoms() as temporary_file:
        rigid = RigidBody.from_pdb(temporary_file.name)
        assert len(rigid) == 10


def test_copy():
    rigid = generate_rigid_from_pdb()
    rigid_copy = rigid.copy()
    assert isinstance(rigid_copy, RigidBody)
    assert len(rigid_copy) == len(rigid)
    assert rigid_copy.coordinates.shape == rigid.coordinates.shape


def test_copy_does_not_contain_any_reference():
    rigid = generate_rigid_from_pdb()
    rigid_copy = rigid.copy()
    rigid.coordinates.fill(0)
    assert rigid.coordinates.all() == 0
    assert not (rigid_copy.coordinates == 0).all()
