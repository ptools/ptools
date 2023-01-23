"""test_spelling - Tests for ptools.spelling."""

from ptools import spelling


def test_plural_simple():
    assert spelling.pluralize("name") == "names"
    assert spelling.pluralize("element") == "elements"


def test_plural_irregular():
    assert spelling.pluralize("index") == "indices"


def test_plural_multiple_words():
    assert spelling.pluralize("residue_name") == "residue_names"
    assert spelling.pluralize("residue_index") == "residue_indices"


def test_plural_single_letter():
    assert spelling.pluralize("x") == "xs"
