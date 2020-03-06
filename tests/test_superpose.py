
import random
import unittest

from pyptools.rigidbody import RigidBody
import pyptools.superpose as superpose


from . import TEST_LIGAND


class TestSuperpose(unittest.TestCase):
    def setUp(self):
        self.reference = RigidBody(TEST_LIGAND)

