
import unittest

from pyptools.rigidbody import AttractRigidBody
from pyptools.forcefield import ForceField, AttractForceField1


from .pyattract import TEST_AMINON, TEST_LIGAND_RED, TEST_RECEPTOR_RED



class TestForceField(unittest.TestCase):

    def test_energy_is_zero(self):
        ff = ForceField(receptor=None, ligand=None)
        self.assertEqual(ff.energy(), 0.0)


class TestAttractForceField1(unittest.TestCase):

    def setUp(self):
        self.receptor = AttractRigidBody(TEST_RECEPTOR_RED)
        self.ligand = AttractRigidBody(TEST_LIGAND_RED)
        self.ff = AttractForceField1(self.receptor, self.ligand, TEST_AMINON)

    def test_foo(self):
        self.ff.energy()
