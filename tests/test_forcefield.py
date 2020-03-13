
import unittest

from pyptools.rigidbody import AttractRigidBody
from pyptools.forcefield import ForceField, AttractForceField1

from .pyattract import TEST_AMINON, TEST_LIGAND_RED, TEST_RECEPTOR_RED
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
        self.ff = AttractForceField1(self.receptor, self.ligand)

    def test_read_ff_params(self):
        ff = AttractForceField1(self.receptor, self.ligand,
                                paramfile=TEST_AMINON)
        assert_array_equal(ff._repulsive, self.ff._repulsive)
        assert_array_equal(ff._attractive, self.ff._attractive)

    def test_calculate_energy(self):
        self.ff.update()
        # Reference values calculated from PTools 1d4b930.
        self.assertAlmostEqual(self.ff.vdw_energy(), -4.85626395114)
        self.assertAlmostEqual(self.ff.electrostatic_energy(), 0.0)


    def test_calculate_energy_with_electrostatic(self):
        self.receptor.atom_charges = np.ones(self.receptor.size())
        self.ligand.atom_charges = np.ones(self.receptor.size()) * 2
        self.ff.update()
        self.assertAlmostEqual(self.ff.vdw_energy(), -4.85626395114)
        self.assertNotEqual(self.ff.electrostatic_energy(), 0.0)


