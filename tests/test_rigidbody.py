"""test_rigidbody - Tests for `ptools.rigidbody module."""

import unittest

from ptools.rigidbody import RigidBody, AttractRigidBody
from ptools.io.pdb import InvalidPDBFormatError

from .testing.io import mk_pdb_10_atoms
from .testing.io.red import mk_red, mk_red_invalid_charges, TestREDBuilder
from .testing.moreassert import assert_array_equal

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


# class TestRigidBody(unittest.TestCase):
#     def setUp(self):
#         with mk_pdb_10_atoms() as tmp_pdb:
#             self.rigid = RigidBody.from_pdb(tmp_pdb.name)

#     def test_constructor(self):
#         self.assertEqual(len(self.rigid), 10)

#     def _assert_copy_successful(self, thecopy):
#         # Check that both RigidBody instances have the same number of atoms.
#         self.assertEqual(len(thecopy), len(self.rigid))
#         self.assertEqual(thecopy.coordinates.shape, (10, 3))

#         # Change parent ParticleCollection coordinates and make sure it does not
#         # affect the copy.
#         ref_coords = self.rigid.coordinates.copy()
#         self.rigid.coordinates.fill(0)
#         assert_array_equal(thecopy.coordinates, ref_coords)

#     def test_copy_constructor1(self):
#         thecopy = self.rigid.copy()
#         self._assert_copy_successful(thecopy)

#     def test_copy_constructor2(self):
#         thecopy = RigidBody.copy(self.rigid)
#         self._assert_copy_successful(thecopy)

#     def test_getitem_returns_rigidbody(self):
#         atoms = self.rigid[:5]
#         self.assertIsInstance(atoms, RigidBody)


# class TestAttractRigidBody(unittest.TestCase):
#     def test_constructor(self):
#         with mk_red() as tmpfile:
#             arb = AttractRigidBody.from_pdb(tmpfile.name)

#         self.assertEqual(len(arb), 10)
#         # !! atom categories stored in AttractRigidBody are array indices,
#         # !! therefore minored by one compared to what is in the RED file.
#         assert_array_equal(arb.atom_categories + 1, TestREDBuilder.categories())
#         assert_array_equal(arb.atom_charges, TestREDBuilder.charges())

#     def test_constructor_fails_no_categories(self):
#         with mk_red(has_categories=False) as tmpfile:
#             err = "Expected atom categories and charges, found"
#             with self.assertRaisesRegex(InvalidPDBFormatError, err):
#                 AttractRigidBody.from_pdb(tmpfile.name)

#     def test_constructor_fails_no_charges(self):
#         with mk_red(has_charges=False) as tmpfile:
#             err = "Expected atom categories and charges, found"
#             with self.assertRaisesRegex(InvalidPDBFormatError, err):
#                 AttractRigidBody.from_pdb(tmpfile.name)

#     def test_constructor_fails_invalid_charges(self):
#         with mk_red_invalid_charges() as tmpfile:
#             err = "Atom charge expects a float"
#             with self.assertRaisesRegex(InvalidPDBFormatError, err):
#                 AttractRigidBody.from_pdb(tmpfile.name)

#     # Ignores R0201: Method could be a function (no-self-use)
#     # pylint: disable=R0201
#     def test_reset_forces(self):
#         with mk_red() as tmpfile:
#             arb = AttractRigidBody.from_pdb(tmpfile.name)
#         assert_array_equal(arb.atom_forces, 0.0)
#         arb.atom_forces[0] = [1, 2, 3]
#         assert_array_equal(arb.atom_forces[0], [1, 2, 3])
#         arb.reset_forces()
#         assert_array_equal(arb.atom_forces, 0.0)
