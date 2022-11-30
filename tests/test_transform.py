"""test_transform.py - Tests for ``ptools.transform``."""

import numpy as np

from .testing import assert_array_almost_equal, assert_array_not_almost_equal
from .testing.dummy import generate_dummy_atomcollection

from ptools import measure, transform

def test_translate():
    atoms = generate_dummy_atomcollection()
    origin = (0, 0, 0)
    center = measure.centroid(atoms)
    transform.translate(atoms, origin - center)
    assert_array_almost_equal(measure.centroid(atoms), origin)


def test_translate_scalar():
    atoms = generate_dummy_atomcollection()
    centroid = measure.centroid(atoms)
    assert_array_almost_equal(centroid, centroid[0])
    scalar = -centroid[0]
    transform.translate(atoms, scalar)
    assert_array_almost_equal(measure.centroid(atoms), (0, 0, 0))


def test_center_without_weigths():
    atoms = generate_dummy_atomcollection()
    origin = np.zeros(3)
    for origin in (np.zeros(3), np.ones(3)):
        assert_array_not_almost_equal(measure.centroid(atoms), origin)
        transform.center_to_origin(atoms, origin)
        assert_array_almost_equal(measure.centroid(atoms), origin)
