"""test_attractrigidbody - Tests for `ptools.rigidbody.AttractRigidBody."""

import numpy as np
import pytest
from pytest import approx

from ptools import AttractRigidBody
from ptools.io import read_attract_topology
from ptools.io.exceptions import InvalidREDFormatError

from .generators import generate_red_file
from .generators.red import RedFileBuilder

# == Tests for AttractRigidBody initialization =========================================


def test_initialization_from_pdb():
    with generate_red_file() as temporary_file:
        rigid = read_attract_topology(temporary_file.name)

    assert isinstance(rigid, AttractRigidBody)
    assert len(rigid) == 10
    assert hasattr(rigid, "typeids")
    assert hasattr(rigid, "charges")
    assert hasattr(rigid, "forces")

    assert rigid.typeids == approx(np.array(RedFileBuilder.typeids()) - 1)
    assert rigid.charges == approx(RedFileBuilder.charges())
    assert rigid.radii == approx(RedFileBuilder.radii())
    assert rigid.forces == approx(np.zeros((10, 3)))


def test_initialization_from_pdb_fails_no_typeids():
    with generate_red_file(has_typeids=False) as temporary_file:
        err = "Expected atom type ids and charges, found"
        with pytest.raises(InvalidREDFormatError, match=err):
            read_attract_topology(temporary_file.name)


def test_initialization_from_pdb_fails_no_charges():
    with generate_red_file(has_charges=False) as temporary_file:
        err = "Expected atom type ids and charges, found"
        with pytest.raises(InvalidREDFormatError, match=err):
            read_attract_topology(temporary_file.name)


def test_constructor_fails_invalid_charges():
    with generate_red_file(invalid_charges=True) as temporary_file:
        err = "Atom charge expects a float"
        with pytest.raises(InvalidREDFormatError, match=err):
            read_attract_topology(temporary_file.name)


# =======================================================================================


def test_apply_forces():
    with generate_red_file() as temporary_file:
        rigid = read_attract_topology(temporary_file.name)

    # Forces should be initialized to zero.
    expected = np.zeros((rigid.size(), 3))
    assert rigid.forces == approx(expected)

    # Applying forces should update the forces array.
    expected = np.ones((rigid.size(), 3))
    rigid.apply_forces(expected)
    assert rigid.forces == approx(expected)

    # Applying forces one more time should add the forces.
    rigid.apply_forces(expected)
    assert rigid.forces == approx(expected * 2)


def test_reset_forces():
    with generate_red_file() as temporary_file:
        rigid = read_attract_topology(temporary_file.name)

    assert rigid.forces == approx(0.0)

    rigid.forces[0] = [1, 2, 3]
    assert rigid.forces[0] == approx([1, 2, 3])

    rigid.reset_forces()
    assert rigid.forces == approx(0.0)
