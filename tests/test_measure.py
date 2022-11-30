"""test_measure.py - Tests for ptools.measure."""

import numpy as np

import ptools
from ptools import measure

from pytest import approx

from .testing.dummy import generate_dummy_atomcollection, generate_balloon
from .testing import assert_array_almost_equal


def test_distance_to_axis():
    obj = generate_balloon((0, 0, 0))

    axis = np.array((1, 0, 0))
    assert measure.distance_to_axis(obj, axis) == approx(0.0)

    obj.coordinates = (0, 1, 0)
    assert measure.distance_to_axis(obj, axis) == approx(1.0)

    obj.coordinates = (0, 1, 1)
    assert measure.distance_to_axis(obj, axis) == approx(2.0 ** 0.5)


def test_centroid():
    atoms = generate_dummy_atomcollection()
    center = measure.centroid(atoms)
    assert_array_almost_equal(center, [4.5, 4.5, 4.5])


def test_center_of_masses():
    # Masses are 1. COM should be same as center
    atoms = generate_dummy_atomcollection()
    center = measure.center_of_mass(atoms)
    assert_array_almost_equal(center, [4.5, 4.5, 4.5])

    # Changes atoms name, therefore element, therefore mass.
    for i in range(5):
        atoms[i].name = "CA"
    for i in range(5, 10):
        atoms[i].name = "NZ"
    atoms.guess_masses()

    center = measure.center_of_mass(atoms)
    assert_array_almost_equal(center, [4.69179, 4.69179, 4.69179])


def test_inertia_tensor():
    from . import TEST_LIGAND
    atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)

    # Reference calculated with MDAnalysis 0.20.1:
    # >>> MDAnalysis.Universe("ligand.pdb").select_atoms("all").moment_of_inertia()
    ref = [
        [3679339.47775172, 694837.16289436, -263651.10452372],
        [694837.16289436, 3803047.59374612, -194611.71739629],
        [-263651.10452372, -194611.71739629, 3425042.27240564],
    ]
    I = measure.inertia_tensor(atoms)
    assert_array_almost_equal(I, ref, decimal=2)


def test_principal_axes():
    from . import TEST_LIGAND
    atoms = ptools.io.pdb.read_pdb(TEST_LIGAND)

    # Calculated with MDAnalysis 0.20.1:
    # >>> MDAnalysis.Universe("ligand.pdb").select_atoms("all").principal_axes()
    ref = [
        [0.65682984, 0.70033642, -0.27946997],
        [0.04052252, 0.33731064, 0.94052084],
        [0.7529492, -0.62908698, 0.1931763],
    ]
    I = measure.principal_axes(atoms)
    assert_array_almost_equal(I, ref)


def test_radius_of_gyration():
    atoms = generate_dummy_atomcollection()
    rgyr = measure.radius_of_gyration(atoms)
    # Reference value calculated with VMD.
    assert rgyr == approx(4.9749369621276855, 1e-6)