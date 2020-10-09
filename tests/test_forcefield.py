
import unittest

from ptools.rigidbody import AttractRigidBody
from ptools.forcefield import ForceField, AttractForceField1

from .attract import TEST_AMINON, TEST_LIGAND_RED, TEST_RECEPTOR_RED
from .testing.moreassert import assert_array_equal

import numpy as np


class TestForceField(unittest.TestCase):

    def test_energy_is_zero(self):
        ff = ForceField(receptor=None, ligand=None)
        self.assertEqual(ff.energy(), 0.0)


class TestAttractForceField1(unittest.TestCase):

    def setUp(self):
        self.receptor = AttractRigidBody(TEST_RECEPTOR_RED)
        self.ligand = AttractRigidBody(TEST_LIGAND_RED)
        self.ff = AttractForceField1(self.receptor, self.ligand, cutoff=5.0)

    def test_read_ff_params(self):
        ff = AttractForceField1(self.receptor, self.ligand,
                                paramfile=TEST_AMINON, cutoff=5.0)
        assert_array_equal(ff._repulsive_parameters, self.ff._repulsive_parameters)
        assert_array_equal(ff._attractive_parameters, self.ff._attractive_parameters)

    def test_calculate_energy(self):
        # Reference values calculated from PTools 1d4b930.
        self.ff.update()
        self.assertAlmostEqual(self.ff.vdw_energy(), -4.85626395114)
        self.assertAlmostEqual(self.ff.electrostatic_energy(), 0.0)

    def test_calculate_energy_with_electrostatic(self):
        # Reference values calculated from PTools 8439c40.
        self.receptor.atom_charges = np.ones(self.receptor.size())
        self.ligand.atom_charges = np.ones(self.ligand.size()) * 2
        self.ff.update()
        self.assertAlmostEqual(self.ff.vdw_energy(), -4.85626395114)
        self.assertAlmostEqual(self.ff.electrostatic_energy(), 77.90251420616622)

    def test_non_bonded_energy_small_cutoff(self):
        # Reference values calculated from PTools 8439c40.
        self.ff.cutoff = 5.0
        self.assertAlmostEqual(self.ff.non_bonded_energy(), -4.85626395114)

    def test_non_bonded_energy_large_cutoff(self):
        # Reference values calculated from PTools 8439c40.
        self.ff.cutoff = 50.0
        self.assertAlmostEqual(self.ff.non_bonded_energy(), -56.10729372698934)

    def test_non_bonded_energy_methods_return_same_result(self):
        # Reference values calculated from PTools 8439c40.
        self.assertAlmostEqual(self.ff.non_bonded_energy(), -4.85626395114)
        self.assertAlmostEqual(self.ff._AttractForceField1__nb_energy_small_cutoff(), -4.85626395114)
        self.assertAlmostEqual(self.ff._AttractForceField1__nb_energy_large_cutoff(), -4.85626395114)
