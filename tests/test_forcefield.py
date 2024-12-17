"""Test for ptools.forcefield module."""

import unittest
from pathlib import Path

import numpy as np

from ptools.attract import AttractRigidBody
from ptools.forcefield import AttractForceField1
from ptools.io import read_attract_topology

from .attract import TEST_LIGAND_RED, TEST_RECEPTOR_RED
from .testing.moreassert import assert_array_almost_equal


# Ignores W0212: Access to a protected member of a client class
# pylint: disable=W0212
class TestAttractForceField1DummyRigid(unittest.TestCase):
    """Tests for AttractForceField1 that do not required an actual AttractRigidBody."""

    path_ff_parameters = Path(__file__).parent / "data" / "ff_parameters.np"

    def setUp(self):
        self.receptor = read_attract_topology(TEST_RECEPTOR_RED)
        self.ligand = read_attract_topology(TEST_LIGAND_RED)
        self.ff = AttractForceField1(self.receptor, self.ligand, cutoff=5.0)

        # Save forcefield parameters for later use in non-regression tests.
        # Uncomment to save new file.
        # self._save_forcefield_parameters()

    def _save_forcefield_parameters(self):
        """Save forcefield parameters for later use in regression tests."""
        with self.path_ff_parameters.open("wb") as f:
            np.save(f, self.ff._attractive_parameters)
            np.save(f, self.ff._repulsive_parameters)

    def _load_forcefield_parameters(self) -> tuple[np.ndarray, np.ndarray]:
        """Load forcefield parameters for regression tests."""
        with self.path_ff_parameters.open("rb") as f:
            attractive = np.load(f)
            repulsive = np.load(f)
        return (attractive, repulsive)

    def test_regression_ff_parameters(self):
        """Regression test for forcefield parameters."""
        attractive, repulsive = self._load_forcefield_parameters()
        assert_array_almost_equal(self.ff._attractive_parameters, attractive)
        assert_array_almost_equal(self.ff._repulsive_parameters, repulsive)


class TestAttractForceField1(unittest.TestCase):
    """Tests for AttractForceField that requires actual AttractRigidBody instances."""

    def setUp(self):
        self.receptor = read_attract_topology(TEST_RECEPTOR_RED)
        self.ligand = read_attract_topology(TEST_LIGAND_RED)
        self.ff = AttractForceField1(self.receptor, self.ligand, cutoff=5.0)

    def test_calculate_energy(self):
        # Reference values calculated from PTools 1d4b930.
        self.ff.update()
        self.assertAlmostEqual(self.ff.vdw_energy(), -4.85626395114)
        self.assertAlmostEqual(self.ff.electrostatic_energy(), 0.0)

    def test_calculate_energy_with_electrostatic(self):
        # Reference values calculated from PTools 8439c40.
        self.receptor.charges = np.ones(self.receptor.size())
        self.ligand.charges = np.ones(self.ligand.size()) * 2
        self.ff.update()
        self.assertAlmostEqual(self.ff.vdw_energy(), -4.85626395114)
        self.assertAlmostEqual(self.ff.electrostatic_energy(), 77.90251420616622)

    def test_non_bonded_energy_small_cutoff(self):
        # Reference values calculated from PTools 8439c40.
        self.ff.cutoff = 5.0
        print(f"{self.ff.non_bonded_energy()=}")
        self.assertAlmostEqual(self.ff.non_bonded_energy(), -4.85626395114)

    def test_non_bonded_energy_large_cutoff(self):
        # Reference values calculated from PTools 8439c40.
        self.ff.cutoff = 50.0
        self.assertAlmostEqual(self.ff.non_bonded_energy(), -56.10729372698934)

    # Ignores W0212: Access to a protected member of a client class
    # pylint: disable=W0212
    def test_non_bonded_energy_methods_return_same_result(self):
        # Reference values calculated from PTools 8439c40.
        self.assertAlmostEqual(self.ff.non_bonded_energy(), -4.85626395114)
        self.assertAlmostEqual(
            self.ff._AttractForceField1__nb_energy_small_cutoff(), -4.85626395114
        )
        self.assertAlmostEqual(
            self.ff._AttractForceField1__nb_energy_large_cutoff(), -4.85626395114
        )
