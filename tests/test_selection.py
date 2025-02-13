import unittest
from pathlib import Path

import numpy as np

from ptools.io import read_pdb

PDB_TEST_SELECTION = Path(__file__).parent / "data" / "test_selection.pdb"


# class TestSelectionSyntax(unittest.TestCase):
#
#     def test_unknown_token_error(self):
#         token = "not_a_valid_token"
#         selection_str = f"{token} A B C"
#         with self.assertRaises(UnknownTokenError) as context:
#             select(selection_str)
#             assert context.exception.token == token


class TestSelectionBase(unittest.TestCase):
    """Base class for selection test cases."""

    def setUp(self):
        self.atoms = read_pdb(PDB_TEST_SELECTION)
        assert len(self.atoms) == 66


class TestSelectionAtomIndex(TestSelectionBase):
    def test_selection_atom_index(self):
        result = self.atoms.select("index 1")
        assert result == self.atoms[:1]

    def test_selection_atom_multiple_indices(self):
        result = self.atoms.select("index 1 2 3")
        assert result == self.atoms[:3]

    def test_selection_atom_range(self):
        result = self.atoms.select("index 1 to 3")
        assert result == self.atoms[:3]

        result = self.atoms.select("index 2:5")
        assert result == self.atoms[1:5]


class TestSelectionResidueIndex(TestSelectionBase):
    def test_selection_residue_index(self):
        result = self.atoms.select("resid 2")
        assert len(result) == 3
        for atom in result:
            assert atom.residue_index == 2

    def test_selection_residue_multiple_indices(self):
        result = self.atoms.select("resid 2 3")
        assert len(result) == 6
        for atom in result:
            assert atom.residue_index in (2, 3)

    def test_selection_residue_range(self):
        result = self.atoms.select("resid 2 to 3")
        assert len(result) == 6
        for atom in result:
            assert atom.residue_index in (2, 3)

        result = self.atoms.select("resid 2:3")
        assert len(result) == 6
        for atom in result:
            assert atom.residue_index in (2, 3)


class TestSelectionAtomName(TestSelectionBase):
    def test_selection_single_atom_name(self):
        result = self.atoms.select("name CA")
        assert len(result) == 29
        for atom in result:
            assert atom.name == "CA"

    def test_selection_multiple_atom_names(self):
        result = self.atoms.select("name CA CB")
        assert len(result) == 36
        for atom in result:
            assert atom.name in ("CA", "CB")


class TestSelectionResidueName(TestSelectionBase):
    def test_selection_single_residue_name(self):
        result = self.atoms.select("resname ALA")
        assert len(result) == 4
        for atom in result:
            assert atom.residue_name == "ALA"

    def test_selection_multiple_residue_names(self):
        result = self.atoms.select("resname ALA ARG GLU")
        assert len(result) == 16
        for atom in result:
            assert atom.residue_name in ("ALA", "ARG", "GLU")


class TestSelectionChainName(TestSelectionBase):
    def test_selection_single_chain(self):
        result = self.atoms.select("chain A")
        assert len(result) == 21
        for atom in result:
            assert atom.chain == "A"

    def test_selection_multiple_chains(self):
        result = self.atoms.select("chain A B")
        self.assertEqual(len(result), 43)
        for atom in result:
            assert atom.chain in ("A", "B")


# == Operator Tests ===============================================================


class TestSelectionAndOperator(TestSelectionBase):
    def test_selection_and_simple(self):
        atoms = self.atoms.select("resid 1 to 14")
        assert len(atoms) == 29

        atoms = self.atoms.select("resid 1 to 14 and chain A")
        assert len(atoms) == 21

        atoms = self.atoms.select("resid 1 to 14 and chain B")
        assert len(atoms) == 8

        atoms = self.atoms.select("resid 1 to 14 and chain B and name CA")
        assert len(atoms) == 4


class TestSelectionOrOperator(TestSelectionBase):
    def test_selection_or_operator(self):
        result = self.atoms.select("resid 2 or resid 3")
        assert len(result) == 6
        for atom in result:
            assert atom.residue_index in (2, 3)


class TestSelectionNotOperator(TestSelectionBase):
    def test_selection_not_chain(self):
        result = self.atoms.select("not chain B")
        assert len(result) == 44

        result = self.atoms.select("not resname ALA ARG GLU")
        assert len(result) == 50

        result = self.atoms.select("not resname ALA and chain A")
        assert len(result) == 19


class TestSelectionArithmeticOperator(TestSelectionBase):
    def test_selection_arithmetic_operator_syntax(self):
        self.atoms.select("resid < 5")
        self.atoms.select("resid< 5")
        self.atoms.select("resid <5")
        self.atoms.select("resid <5")
        self.atoms.select("resid <= 5")

    def test_selection_less_than(self):
        result = self.atoms.select("resid < 5")
        assert len(result) == 8
        for atom in result:
            assert atom.residue_index < 5

    def test_selection_less_than_or_equal(self):
        result = self.atoms.select("resid <= 5")
        assert len(result) == 10
        for atom in result:
            assert atom.residue_index <= 5


class TestSelectionBooleanProperty(TestSelectionBase):
    def test_selection_boolean_property(self):
        result = self.atoms.select("hetero")
        assert len(result) == 3

        result = self.atoms.select("not hetero")
        assert len(result) == 63


class TestSelectionDynamicProperty(TestSelectionBase):
    def test_selection_dynamic_property(self):
        self.atoms.add_atom_property(
            "cherry", "cherries", np.full(len(self.atoms), "foo")
        )

        result = self.atoms.select("cherry foo")
        assert len(result) == 66


class TestSelectionKeywords(TestSelectionBase):
    def test_selection_water(self):
        result = self.atoms.select("water")
        assert len(result) == 0


class TestSelectionByProperty(TestSelectionBase):

    def test_single_value(self):
        result = self.atoms.select_by_property("index", 1)
        assert result == self.atoms[:1]

        result = self.atoms.select_by_property("chain", "A")
        assert len(result) == 21

        result = self.atoms.select_by_property("resid", 2)
        assert len(result) == 3

    def test_list(self):
        result = self.atoms.select_by_property("index", [1, 2, 3])
        assert result == self.atoms[:3]

        result = self.atoms.select_by_property("resid", [5, 6, 7])
        assert len(result) == 7

