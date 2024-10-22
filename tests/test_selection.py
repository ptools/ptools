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
            assert context.exception.token == token


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
        assert len(atoms) == 3
        for atom in atoms:
            assert atom.residue_index == 2

    def test_selection_residue_multiple_indices(self):
        atoms = select("resid 2 3", self.atoms)
        assert len(atoms) == 6
        for atom in atoms:
            assert atom.residue_index in (2, 3)

    def test_selection_residue_range(self):
        atoms = select("resid 2 to 3", self.atoms)
        assert len(atoms) == 6
        for atom in atoms:
            assert atom.residue_index in (2, 3)

        atoms = select("resid 2:3", self.atoms)
        assert len(atoms) == 6
        for atom in atoms:
            assert atom.residue_index in (2, 3)


class TestSelectionAtomName(TestSelectionBase):

    def test_selection_single_atom_name(self):
        atoms = select("name CA", self.atoms)
        assert len(atoms) == 29
        for atom in atoms:
            assert atom.name == "CA"

    def test_selection_multiple_atom_names(self):
        atoms = select("name CA CB", self.atoms)
        assert len(atoms) == 36
        for atom in atoms:
            assert atom.name in ("CA", "CB")


class TestSelectionResidueName(TestSelectionBase):

    def test_selection_single_residue_name(self):
        atoms = select("resname ALA", self.atoms)
        assert len(atoms) == 4
        for atom in atoms:
            assert atom.residue_name == "ALA"

    def test_selection_multiple_residue_names(self):
        atoms = select("resname ALA ARG GLU", self.atoms)
        assert len(atoms) == 16
        for atom in atoms:
            assert atom.residue_name in ("ALA", "ARG", "GLU")


class TestSelectionChainName(TestSelectionBase):

    def test_selection_single_chain(self):
        atoms = select("chain A", self.atoms)
        assert len(atoms) == 21
        for atom in atoms:
            assert atom.chain == "A"

    def test_selection_multiple_chains(self):
        atoms = select("chain A B", self.atoms)
        self.assertEqual(len(atoms), 43)
        for atom in atoms:
            assert atom.chain in ("A", "B")


# == Operator Tests ===============================================================

class TestSelectionAndOperator(TestSelectionBase):

    def test_selection_and_simple(self):
        atoms = select("resid 1 to 14", self.atoms)
        assert len(atoms) == 29

        atoms = select("resid 1 to 14 and chain A", self.atoms)
        assert len(atoms) == 21

        atoms = select("resid 1 to 14 and chain B", self.atoms)
        assert len(atoms) == 8

        atoms = select("resid 1 to 14 and chain B and name CA", self.atoms)
        assert len(atoms) == 4


class TestSelectionOrOperator(TestSelectionBase):

    def test_selection_or_operator(self):
        atoms = select("resid 2 or resid 3", self.atoms)
        assert len(atoms) == 6
        for atom in atoms:
            assert atom.residue_index in (2, 3)


class TestSelectionNotOperator(TestSelectionBase):

    def test_selection_not_chain(self):
        atoms = select("not chain B", self.atoms)
        assert len(atoms) == 44

        atoms = select("not resname ALA ARG GLU", self.atoms)
        assert len(atoms) == 50

        atoms = select("not resname ALA and chain A", self.atoms)
        assert len(atoms) == 19

