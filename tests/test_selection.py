
import unittest
from pathlib import Path

import pytest

from ptools.io import read_pdb
from ptools.selection import MisformattedExpressionError, select, UnknownTokenError

PDB_TEST_SELECTION = Path(__file__).parent / "data" / "test_selection.pdb"


class TestSelectionSyntax(unittest.TestCase):

    def test_unknown_token_error(self):
        token = "not_a_valid_token"
        selection_str = f"{token} A B C"
        with self.assertRaises(UnknownTokenError) as context:
            select(selection_str)
            self.assertEqual(context.exception.token, token)


class TestSelectionBase(unittest.TestCase):
    """Base class for selection test cases."""

    def setUp(self):
        self.atoms = read_pdb(PDB_TEST_SELECTION)
        assert len(self.atoms) == 66


class TestSelectionAtomIndex(TestSelectionBase):

    def test_selection_atom_index(self):
        atoms = select("index 1", self.atoms)
        self.assertEqual(len(atoms), 1)
        self.assertEqual(atoms[0].index, 1)

    def test_selection_atom_multiple_indices(self):
        atoms = select("index 1 2 3", self.atoms)
        assert atoms == self.atoms[:3]

    def test_selection_atom_range(self):
        atoms = select("index 1 to 3", self.atoms)
        assert atoms == self.atoms[:3]

        atoms = select("index 2:5", self.atoms)
        assert atoms == self.atoms[1:5]


class TestSelectionResidueIndex(TestSelectionBase):

    def test_selection_residue_index(self):
        atoms = select("resid 2", self.atoms)
        self.assertEqual(len(atoms), 3)
        for atom in atoms:
            self.assertEqual(atom.residue_index, 2)

    def test_selection_residue_multiple_indices(self):
        atoms = select("resid 2 3", self.atoms)
        self.assertEqual(len(atoms), 6)
        for atom in atoms:
            self.assertIn(atom.residue_index, (2, 3))

    def test_selection_residue_range(self):
        atoms = select("resid 2 to 3", self.atoms)
        self.assertEqual(len(atoms), 6)
        for atom in atoms:
            self.assertIn(atom.residue_index, (2, 3))

        atoms = select("resid 2:3", self.atoms)
        self.assertEqual(len(atoms), 6)
        for atom in atoms:
            self.assertIn(atom.residue_index, (2, 3))


# class TestSelectionRange(TestSelectionBase):
#
#     def test_selection_range(self):
#         atoms = select("index 1 to 5", self.atoms)
#         self.assertEqual(len(atoms), 5)
#         for i, atom in enumerate(atoms):
#             self.assertEqual(atom.index, i + 1)
#


# == DEPRECATED ===================================================================

# Parametrized tests for selection by atom attribute.
#
# Tuples are in the form (<selection token>, <values>, <expected number of particles>).

# @pytest.mark.parametrize(
#     "attr, values, n_atoms_expected",
#     [
#         ("index", (1, ), 1),
#         ("index", (1, 2), 2),
#         ("resid", (2, ), 3),
#         ("resid", (2, 3), 6),
#         ("name", ("CA", ), 29),
#         ("name", ("CA", "CB"), 36),
#         ("resname", ("ALA", ), 4),
#         ("resname", ("ALA", "ARG", "GLU"), 16),
#         ("chain", ("A", ), 21),
#         ("chain", ("A", "B"), 43),
#     ]
# )
# def test_selection_atom_attr(attr, values, n_atoms_expected):
#     """Generic function for testing selection by atom attribute."""
#     atoms = read_pdb(PDB_TEST_SELECTION)
#     assert len(atoms) == 66
#
#     selection_str = f"{attr} {' '.join(str(v) for v in values)}"
#     atoms = select(selection_str, atoms)
#     assert len(atoms) == n_atoms_expected
#     for atom in atoms:
#         assert getattr(atom, attr) in values

# == DEPRECATED - END ============================================================


# == Operator Tests ===============================================================

# class TestSelectionAndOperator(TestSelectionBase):
#
#     def test_selection_and_simple(self):
#         atoms = select("resid 1 to 14", self.atoms)
#         self.assertEqual(len(atoms), 29)
#
#         atoms = select("resid 1 to 14 and chain A", self.atoms)
#         self.assertEqual(len(atoms), 21)
#
#         atoms = select("resid 1 to 14 and chain B", self.atoms)
#         self.assertEqual(len(atoms), 8)
#
#     def test_selection_and_multiple(self):
#         atoms = select("resid 1 to 14 and chain B and name CA", self.atoms)
#         self.assertEqual(len(atoms), 4)
#
#
# class TestSelectionOrOperator(TestSelectionBase):
#
#     def test_selection_or_operator(self):
#         atoms = select("resid 2 or resid 3", self.atoms)
#         self.assertEqual(len(atoms), 6)
#         for atom in atoms:
#             self.assertIn(atom.resid, (2, 3))


# class TestSelectionNotOperator(TestSelectionBase):

#     def test_selection_not_chain(self):
#         atoms = select("not chain B", self.atoms)
#         self.assertEqual(len(atoms), 0)

#     def test_selection_not_resname(self):
#         atoms = select("not resname ALA ARG GLU", self.atoms)
#         self.assertEqual(len(atoms), 772)
