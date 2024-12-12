"""test_atomattrs - Tests for `ptools.atom.AtomAttrs`."""

import pytest
from pytest import approx

from ptools.array3d import Invalid3DArrayError
from ptools.atomattrs import AtomAttrs

from .generators import generate_atom_attrs
from .testing.moreassert import assert_array_almost_equal, assert_array_not_almost_equal


def test_empty_initializer():
    atom = AtomAttrs()
    assert atom.name == "XXX"
    assert atom.index == 0
    assert atom.residue_name == "XXX"
    assert atom.residue_index == 0
    assert atom.chain == "X"
    assert atom.coordinates == approx((0, 0, 0))


def test_initialize_with_bad_coordinates():
    with pytest.raises(Invalid3DArrayError) as excinfo:
        AtomAttrs(coordinates=(1, 2))

    expected = "cannot initialize 3D-coordinates from array with shape (2,)"
    assert str(excinfo.value) == expected


def test_set_bad_coordinates():
    atom = AtomAttrs()
    with pytest.raises(Invalid3DArrayError) as excinfo:
        atom.coordinates = (1, 2)
    expected = "cannot initialize 3D-coordinates from array with shape (2,)"
    assert str(excinfo.value) == expected


def test_set_correct_coordinates():
    atom = AtomAttrs()
    assert atom.coordinates == approx((0, 0, 0))
    atom.coordinates = (1, 2, 3)
    assert atom.coordinates == approx((1, 2, 3))


def test_equality_operator():
    left, right = AtomAttrs(), AtomAttrs()
    assert left == right

    left.name = right.name + "X"
    assert left != right

    left.name = right.name
    left.coordinates = right.coordinates + 1
    assert left != right


def test_copy():
    left = generate_atom_attrs()
    right = left.copy()

    # Checks equality in the standard fashion for all attributes expect ``coordinates``.
    for attr in left.__dict__.keys() - ("coordinates",):
        assert getattr(left, attr) == getattr(right, attr)

    # Checks numpy array of coordinates is actually a copy.
    assert left.coordinates is not right.coordinates
    assert left.coordinates == approx(right.coordinates)

    source_coordinates = left.coordinates.copy()
    left.coordinates += 12

    # right coordinates should remain unchanged.
    assert right.coordinates == approx(source_coordinates)

    # ...while original coordinates should have changed.
    assert_array_not_almost_equal(left.coordinates, source_coordinates)


def test_guess_atom():
    atom = AtomAttrs()
    assert atom.element == "X"  # default value
    atom.name = "CA"
    atom.guess_element()
    assert atom.element == "C"


def test_guess_mass():
    atom = AtomAttrs()
    assert atom.mass == 1.0  # default value
    atom.name = "CA"
    atom.guess_element()
    atom.guess_mass()
    assert atom.mass == approx(12.011)
