"""test_spelling - Tests for ptools.spelling."""

from ptools import spelling


def test_plural_simple():
    assert spelling.plural("name") == "names"
    assert spelling.plural("element") == "elements"


def test_plural_irregular():
    assert spelling.plural("index") == "indices"


def test_plural_multiple_words():
    assert spelling.plural("residue_name") == "residue_names"
    assert spelling.plural("residue_index") == "residue_indices"


def test_plural_single_letter():
    assert spelling.plural("x") == "xs"
